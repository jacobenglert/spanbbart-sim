# Program Name: soft-spanbbart-sim-production.R
# Author:       Jacob Englert
# Date:         27 January 2024
# Purpose:      Simulate overdispersed counts from a complicated mean function
#               and attempt to recover using soft spannbart


# Load Packages -----------------------------------------------------------
library(tidyverse)
library(nbbart)
library(sf)


# Get Parameters ----------------------------------------------------------

# Capture command line arguments
args <- commandArgs(trailingOnly = TRUE)
key <- as.numeric(args[1])
seed <- as.numeric(args[2])
ID <- paste0(sprintf("%03d", key), '-', sprintf("%03d", seed))

# Read parameters from csv file
params <- read_csv(here::here('Params','params-new.csv'), show_col_types = F)

# Subset parameters using key
params <- params[key,-1]

# Set parameters in the current R session
for(i in 1:ncol(params)) {
  assign(names(params[,i]), eval(parse(text = params[,i,drop = TRUE])))
}


# Import Reference Data ---------------------------------------------------

# ATL Zip Code Shapefile (with Population)
data <- read_rds(here::here('Misc','GAZIP05.rds')) |>
  sf::st_as_sf()

# Observed ATL 05-07 predictor correlations
data_cor <- read_rds(here::here('Misc','obs_cor.rds'))


# Simulate Data -----------------------------------------------------------

set.seed(seed)
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
x2 <- apply(x2, 2, \(x) (x - min(x)) / (max(x) - min(x))) #+ 1e-4
colnames(x2) <- paste0('X2_', 1:nx2)

# Fixed parameters
set.seed(seed) # reset seed
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
fit <- soft_spanbbart(x1 = x1, x2 = x2, y = y, s = s, t = t,
                      offset = offset, geo = geo, seed = 1, light_store = TRUE,
                      num_iter = num_iter, num_burn = num_burn, num_thin = num_thin, 
                      m = m, k = k, base = base, power = power,
                      soft = soft, sparse = sparse)


# Compute ALE -------------------------------------------------------------

# Custom BART prediction function (returns entire posterior distribution)
pred_fun <- function(model, newdata) sapply(seq.int(num_burn + 1, num_iter, num_thin), \(i) model$predict_iteration(newdata, i))

# Compute first-order ALE for each predictor
ale1 <- parallel::mclapply(1:nx2, \(j) mcmc_ale(X = x2, model = fit$bart, 
                                                pred_fun = pred_fun, J = j, K = K,
                                                f_true = f, 
                                                center = 'median'), 
                           mc.cores = 5)

ale1 <- do.call(rbind, ale1) |>
  mutate(ID = ID)

# Compute second-order ALEs for first 5 predictors
pairs <- combn(1:5, 2, simplify = FALSE)
ale2 <- parallel::mclapply(pairs, \(j) mcmc_ale(X = x2, model = fit$bart, 
                                                pred_fun = pred_fun, J = j, K = K,
                                                f_true = f, 
                                                center = 'median',
                                                include_main_effects = FALSE),
                           mc.cores = 5)

ale2 <- do.call(rbind, ale2) |>
  mutate(ID = ID)

# Compute second-order + first-order ALEs for first 5 predictors
ale3 <- parallel::mclapply(pairs, \(j) mcmc_ale(X = x2, model = fit$bart, 
                                                pred_fun = pred_fun, J = j, K = K,
                                                f_true = f, 
                                                center = 'median',
                                                include_main_effects = TRUE),
                           mc.cores = 5)

ale3 <- do.call(rbind, ale3) |>
  mutate(ID = ID)


# Compute Simulation Statistics -------------------------------------------

# Implicit intercept
alpha_true <- alpha + mean(f(x2)) + mean(nu) + log(xi)
alpha_est <- median(fit$logmean)# mean(fit$alpha)
alpha_bias <- alpha_est - alpha_true
alpha_lower <- quantile(fit$alpha, 0.025)
alpha_upper <- quantile(fit$alpha, 0.975)
alpha_coverage <- as.numeric(alpha_true >= alpha_lower & alpha_true <= alpha_upper)

# Fixed effects
beta_est <- apply(fit$beta, 2, \(x) quantile(x, 0.500)) #colMeans(fit$beta)
beta_bias <- beta_est - beta
beta_lower <- apply(fit$beta, 2, \(x) quantile(x, 0.025))
beta_upper <- apply(fit$beta, 2, \(x) quantile(x, 0.975))
beta_coverage <- as.numeric(beta >= beta_lower & beta <= beta_upper)

# Dispersion parameter
xi_est <- median(fit$xi)
xi_bias <- xi_est - xi
xi_lower <- quantile(fit$xi, 0.025)
xi_upper <- quantile(fit$xi, 0.975)
xi_coverage <- as.numeric(xi >= xi_lower & xi <= xi_upper)

# Spatial correlation
rho_est <- median(fit$rho)
rho_bias <- rho_est - rho
rho_lower <- quantile(fit$rho, 0.025)
rho_upper <- quantile(fit$rho, 0.975)
rho_coverage <- as.numeric(rho >= rho_lower & rho <= rho_upper)

# Spatial marginal variance
tau2_est <- median(fit$tau2)
tau2_bias <- tau2_est - tau2
tau2_lower <- quantile(fit$tau2, 0.025)
tau2_upper <- quantile(fit$tau2, 0.975)
tau2_coverage <- as.numeric(tau2 >= tau2_lower & tau2 <= tau2_upper)

# Spatial random effects
nu_true <- nu - mean(nu)
nu_bias <- mean(colMeans(fit$nu) - nu_true)
nu_lower <- apply(fit$nu, 2, \(x) quantile(x, 0.025))
nu_upper <- apply(fit$nu, 2, \(x) quantile(x, 0.975))
nu_coverage <- mean(as.numeric(nu_true >= nu_lower & nu_true <= nu_upper))
nu_rmse <- sqrt(mean((colMeans(fit$nu) - nu_true)^2))

# BART component
f_true <- f(x2) + alpha + mean(nu) + log(xi)
bart_mean <- rowMeans(pred_fun(fit$bart, x2))
bart_median <- apply(pred_fun(fit$bart, x2), 1, \(x) quantile(x, 0.500))
bart_lower <- apply(pred_fun(fit$bart, x2), 1, \(x) quantile(x, 0.025))
bart_upper <- apply(pred_fun(fit$bart, x2), 1, \(x) quantile(x, 0.975))
bart_bias <- mean(bart_mean - f_true)
bart_coverage <- mean(as.numeric(f_true >= bart_lower & f_true <= bart_upper))
bart_rmse <- sqrt(mean((bart_mean - f_true)^2))



# Compile Results ---------------------------------------------------------
stats <- data.frame(param = c('alpha',
                              paste0('beta', 1:nx1),
                              'xi',
                              'rho',
                              'tau2',
                              'nu',
                              'G'),
                    est = c(alpha_est, 
                            beta_est, 
                            xi_est, 
                            rho_est, 
                            tau2_est, 
                            NA, 
                            NA),
                    lower = c(alpha_lower, 
                              beta_lower, 
                              xi_lower, 
                              rho_lower, 
                              tau2_lower, 
                              NA, 
                              NA),
                    upper = c(alpha_upper, 
                              beta_upper, 
                              xi_upper, 
                              rho_upper, 
                              tau2_upper, 
                              NA, 
                              NA),
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
  mutate(ID = ID)

results <- list(stats = stats, ale1 = ale1, ale2 = ale2, ale3 = ale3)


# Output Results ----------------------------------------------------------
path <- file.path(here::here(),'..','..','..','projects','hhchang','jrengle','spanbbart-sim','Results','temp')
write_rds(results, paste0(path, '/', ID, '.rds'))
