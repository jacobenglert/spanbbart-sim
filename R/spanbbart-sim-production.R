# Program Name: spanbbart-sim-production.R
# Author:       Jacob Englert
# Date:         28 December 2023
# Purpose:      Simulate overdispersed counts from a complicated mean function
#               and attempt to recover using spannbart


# Load Packages -----------------------------------------------------------
library(tidyverse)
library(nbbart)
library(sf)


# Get Parameters ----------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
names(args) <- c('index')
index <- as.numeric(args['index'])
params <- read_csv(here::here('Params','params.csv'), show_col_types = F)[,-1]
for(i in 1:ncol(params)){
  assign(names(params[index,i]), eval(parse(text = params[index, i, drop = TRUE])))
}


# Import Reference Data ---------------------------------------------------

# ATL Zip Code Shapefile (with Population)
data <- read_rds(here::here('Misc','GAZIP05.rds')) |>
  sf::st_as_sf()

# Observed ATL 05-07 predictor correlations
data_cor <- read_rds(here::here('Misc','obs_cor.rds'))


# Simulate Data -----------------------------------------------------------

set.seed(index)
nT  <- nT             # Length of time-series
nS  <- nrow(data)     # Number of locations
n   <- nT*nS          # Total number of observations

pop <- rep(data$POPULATION, times = nT) # Repeat population vector
geo <- rep(sf::st_geometry(data), times = nT)   # Repeat geometry vector
offset <- log(pop)                      # Create offset

# Confounders
nx1 <- 4
x1 <- t(replicate(n, runif(nx1, 0, 1)))
colnames(x1) <- paste0('X1_', 1:nx1)

# Exposures
nx2 <- 10
nx2a <- 5
nx2b <- 5
# x2 <- t(replicate(n, runif(nx2, 0, 1)))
# x2cov <- diag(10)
x2cov <- rbind(cbind(data_cor, matrix(0, nx2b, nx2b)),
               cbind(matrix(0, nx2a, nx2a), diag(nx2b)))

set.seed(1) # Use the same seed to ensure exposure quantiles do not change
x2 <- mvtnorm::rmvnorm(n, sigma = x2cov)
x2 <- apply(x2, 2, \(x) (x - min(x)) / (max(x) - min(x))) + 1e-4
colnames(x2) <- paste0('X2_', 1:nx2)

# Fixed parameters
set.seed(index) # reset seed
alpha <- -10    # fixed intercept
beta  <- matrix(c(-2, -1, 1, 2), ncol = 1) # fixed effects

# Spatial random effects
# Spatial weight matrices
W <- spdep::poly2nb(data$geometry) |> spdep::nb2mat(style = 'B')
D <- diag(rowSums(W))

tau2  <- 0.3    # marginal variance
rho   <- 0.9    # correlation
nu    <- mvtnorm::rmvnorm(1, sigma = tau2 * solve(D - rho * W))[1,] # effects
x_nu  <- kronecker(matrix(1, nrow = nT, ncol = 1), diag(nS)) # design matrix

# BART component (Friedman 2001)
f <- function(x){
  (10*sin(pi*x[,1]*x[,2]) + 20*(x[,3] - .5)^2 + 10*x[,4] + 5*x[,5]) / 5
}

# Linear predictor
eta <- offset + alpha + x1 %*% beta + f(x2) + x_nu %*% nu

# Dispersion parameter
xi <- 1

# Simulate counts
y <- rnbinom(n, size = xi, prob = 1 / (1 + exp(eta)))


# Fit Model ---------------------------------------------------------------

# Space-time indicators
s <- rep(1:nS, times = nT)
t <- rep(1:nT, each = nS)

# BART model fit
fit <- spanbbart(x1 = x1, x2 = x2, y = y, s = s, t = t,
                 offset = offset, geo = geo, seed = 1, light_store = TRUE,
                 num_iter = num_iter, num_burn = num_burn, num_thin = num_thin, 
                 m = m, k = k, base = base, power = power)


# Compute ALE -------------------------------------------------------------

# Custom BART prediction function (returns entire posterior distribution)
pred_fun <- function(model, newdata) model$predict(newdata)[,c(T,rep(F, num_thin - 1))]

# Compute first-order ALE for each predictor
K <- K
ale1 <- parallel::mclapply(1:nx2, \(j) mcmc_ale(x2, fit$bart, pred_fun, j, K, f_true = f), mc.cores = 5)
names(ale1) <- colnames(x2)

# Create data frame of first-order ALEs
ale1df <- mapply(\(x, n){
  df <- data.frame(var = n, x = x$x, est = x$est, lcl = x$lcl, ucl = x$ucl, truth = x$true)
  return(df)
  }, x = ale1, n = names(ale1), SIMPLIFY = FALSE) |>
  bind_rows() |>
  mutate(ID = index)

# Compute second-order ALEs for first 5 predictors
pairs <- combn(1:5, 2, simplify = FALSE)
ale2 <- parallel::mclapply(pairs, \(j) mcmc_ale(x2, fit$bart, pred_fun, j, K, f_true = f), mc.cores = 5)

# Create data frame of second-order ALEs
ale2df <- mapply(\(x, idx){
  
  n <- colnames(x2)[idx]
  
  x1 <- x$x[[1]]
  x2 <- x$x[[2]]
  
  w1 <- c(diff(x1)[1], diff(x1)) / 2
  w2 <- c(rev(abs(diff(rev(x1)))), rev(abs(diff(x1)))[1]) / 2
  h1 <- c(diff(x2)[1], diff(x2)) / 2
  h2 <- c(rev(abs(diff(rev(x2)))), rev(abs(diff(x2)))[1]) / 2
  
  # Estimate
  f <- as.numeric(x$est)
  truth <- as.numeric(x$true)
  
  f2 <- as.numeric(sweep(sweep(x$est, 1, ale1[[idx[1]]]$est, '+'), 2, ale1[[idx[2]]]$est, '+'))
  truth2 <- as.numeric(sweep(sweep(x$true, 1, ale1[[idx[1]]]$true, '+'), 2, ale1[[idx[2]]]$true, '+'))
  
  df <- data.frame(x1, x2) |>
    expand.grid() |>
    mutate(w1 = rep(w1, times = length(x$x[[2]])),
           w2 = rep(w2, times = length(x$x[[2]])),
           h1 = rep(h1, each = length(x$x[[1]])),
           h2 = rep(h2, each = length(x$x[[1]]))) |>
    mutate(est = f,
           truth = truth,
           lcl = as.numeric(x$lcl),
           ucl = as.numeric(x$ucl),
           est_main = f2,
           truth_main = truth2) |>
    mutate(var1 = n[1], var2 = n[2])
  
  return(df)
}, x = ale2, idx = pairs, SIMPLIFY = FALSE) |>
  bind_rows() |>
  mutate(ID = index)


# Compute Simulation Statistics -------------------------------------------

# Implicit intercept
alpha_true <- alpha + mean(f(x2))
alpha_bias <- mean(fit$alpha) - alpha_true
alpha_lower <- quantile(fit$alpha, 0.025)
alpha_upper <- quantile(fit$alpha, 0.975)
alpha_coverage <- as.numeric(alpha_true >= alpha_lower & alpha_true <= alpha_upper)

# Fixed effects
beta_bias <- colMeans(fit$beta) - beta
beta_lower <- apply(fit$beta, 2, \(x) quantile(x, 0.025))
beta_upper <- apply(fit$beta, 2, \(x) quantile(x, 0.975))
beta_coverage <- as.numeric(beta >= beta_lower & beta <= beta_upper)

# Dispersion parameter
xi_bias <- mean(fit$xi) - xi
xi_lower <- quantile(fit$xi, 0.025)
xi_upper <- quantile(fit$xi, 0.975)
xi_coverage <- as.numeric(xi >= xi_lower & xi <= xi_upper)

# Spatial correlation
rho_bias <- mean(fit$rho) - rho
rho_lower <- quantile(fit$rho, 0.025)
rho_upper <- quantile(fit$rho, 0.975)
rho_coverage <- as.numeric(rho >= rho_lower & rho <= rho_upper)

# Spatial marginal variance
tau2_bias <- mean(fit$tau2) - tau2
tau2_lower <- quantile(fit$tau2, 0.025)
tau2_upper <- quantile(fit$tau2, 0.975)
tau2_coverage <- as.numeric(tau2 >= tau2_lower & tau2 <= tau2_upper)

# Spatial random effects
nu_bias <- mean(colMeans(fit$nu) - nu)
nu_lower <- apply(fit$nu, 2, \(x) quantile(x, 0.025))
nu_upper <- apply(fit$nu, 2, \(x) quantile(x, 0.975))
nu_coverage <- mean(as.numeric(nu >= nu_lower & nu <= nu_upper))
nu_rmse <- sqrt(mean((colMeans(fit$nu) - nu)^2))

# BART component
f_true <- f(x2) + alpha
bart_mean <- rowMeans(fit$bart$predict(x2)[,c(T,rep(F, num_thin - 1))])
bart_median <- apply(fit$bart$predict(x2)[,c(T,rep(F, num_thin - 1))], 1, \(x) quantile(x, 0.500))
bart_lower <- apply(fit$bart$predict(x2)[,c(T,rep(F, num_thin - 1))], 1, \(x) quantile(x, 0.025))
bart_upper <- apply(fit$bart$predict(x2)[,c(T,rep(F, num_thin - 1))], 1, \(x) quantile(x, 0.975))
bart_bias <- mean(bart_mean - f_true)
bart_coverage <- mean(as.numeric(f_true >= bart_lower & f_true <= bart_upper))
bart_rmse <- sqrt(mean((bart_mean - f_true)^2))



# Compile Results ---------------------------------------------------------
sim_results <- data.frame(param = c('alpha',
                                    paste0('beta', 1:nx1),
                                    'xi',
                                    'rho',
                                    'tau2',
                                    'nu',
                                    'G'),
                          bias = c(alpha_bias, 
                                   beta_bias, 
                                   xi_bias, 
                                   rho_bias, 
                                   tau2_bias, 
                                   nu_bias, 
                                   bart_bias),
                          coverage = c(alpha_coverage, 
                                       beta_coverage, 
                                       xi_coverage, 
                                       rho_coverage, 
                                       tau2_coverage, 
                                       nu_coverage, 
                                       bart_coverage),
                          rmse = c(NA,rep(NA, nx1), NA, NA, NA, nu_rmse, bart_rmse)) |>
  mutate(ID = index)

results <- list(stats = sim_results, ale1 = ale1df, ale2 = ale2df)


# Output Results ----------------------------------------------------------
message(paste0("SLURM_ARRAY_TASK_ID is : ", index))
if(!dir.exists(here::here('Results','temp'))) dir.create(here::here('Results','temp'))
write_rds(results, here::here('Results','temp', paste0(sprintf("%04d", index), ".rds")))
