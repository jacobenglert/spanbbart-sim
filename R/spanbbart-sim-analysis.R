# Program Name: spanbbart-sim-analysis.R
# Author:       Jacob Englert
# Date:         2 January 2024
# Purpose:      Analyze results of the spanbbart simulation.


# Load Packages -----------------------------------------------------------
library(tidyverse)


# Load Results ------------------------------------------------------------
sim_date <- '02JAN2024'
results <- read_rds(here::here('Results', paste0(sim_date, '.rds')))


# Summarize Simulation Statistics -----------------------------------------

# Bias, coverage, and RMSE
lower <- function(x) quantile(x, 0.025, na.rm = TRUE)
upper <- function(x) quantile(x, 0.975, na.rm = TRUE)
stats <- lapply(results, '[[', 'stats') |>
  bind_rows() |>
  group_by(param) |>
  summarise(across(c(bias, coverage, rmse), list(mean, lower, upper)))

# First-order ALE
ale1 <- lapply(results, '[[', 'ale1') |>
  bind_rows() |>
  mutate(x_idx = row_number(), .by = ID) |>
  #mutate(x = round(x, .1)) |>
  summarise(across(c(est, lcl, ucl, truth), mean), .by = c(x, var))

ale1 |>
  ggplot(aes(x = x, ymin = lcl, ymax = ucl)) +
  geom_line(aes(y = est)) +
  geom_line(aes(y = truth), col = 'red') +
  geom_ribbon(alpha = 0.4) +
  geom_hline(yintercept = 0, lty = 2) +
  facet_wrap(~factor(var, levels = colnames(x2)), ncol = 5, scales = 'free_y') +
  theme_bw() +
  labs(title = 'ALE Plots of First-Order Main Effects',
       x = 'Predictor Value',
       y = 'ALE')



# Second-order ALE
ale2 <- lapply(results, '[[', 'ale2') |>
  bind_rows() |>
  mutate(x_idx = row_number(), .by = ID) |>
  # mutate(x1 = round(x1, 2),
  #        x2 = round(x2, 2)) |>
  summarise(across(c(est, lcl, ucl, truth, est_main, truth_main), mean), 
            .by = c(x1, x2, w1, w2, h1, h2, var1, var2)) |>
  arrange(var1, var2, x1, x2)

# ale2df <- ale2 |>
#   group_by(var1, var2) |>
#   group_split() |>
#   lapply(\(x){
#     x1 <- unique(x$x1)
#     x2 <- unique(x$x2)
#     w1 <- c(diff(x1)[1], diff(x1)) / 2
#     w2 <- c(rev(abs(diff(rev(x1)))), rev(abs(diff(x1)))[1]) / 2
#     h1 <- c(diff(x2)[1], diff(x2)) / 2
#     h2 <- c(rev(abs(diff(rev(x2)))), rev(abs(diff(x2)))[1]) / 2
#     
#     x |>
#       mutate(w1 = rep(w1, each = length(unique(x$x2))),
#              w2 = rep(w2, each = length(unique(x$x2))),
#              h1 = rep(h1, times = length(unique(x$x1))),
#              h2 = rep(h2, times = length(unique(x$x1))))
#     
#   }) |>
#   bind_rows()


# Estimate (without Main Effects)
limits <- c(min(ale2$truth, ale2$est), max(ale2$truth, ale2$est))
ale2_est <- ale2 |>
  ggplot(aes(xmin = x1 - w1, xmax = x1 + w2, ymin = x2 - h1, ymax = x2 + h2)) +
  geom_rect(aes(fill = est)) +
  theme_bw() +
  coord_fixed() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_gradient2(low = 'blue', high = 'red', limits = limits) +
  facet_grid(var2 ~ var1) +
  theme_bw() +
  labs(title = 'Estimated Second-Order ALEs (excluding Main Effects)',
       x = 'Predictor 1 Value',
       y = 'Predictor 2 Value',
       fill = 'ALE')

# Truth (without Main Effects)
ale2_true <- ale2 |>
  ggplot(aes(xmin = x1 - w1, xmax = x1 + w2, ymin = x2 - h1, ymax = x2 + h2)) +
  geom_rect(aes(fill = truth)) +
  theme_bw() +
  coord_fixed() +
  scale_x_continuous(expand = c(0,0), breaks = seq(0,1,0.5)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_gradient2(low = 'blue', high = 'red', limits = limits) +
  facet_grid(var2 ~ var1) +
  theme_bw() +
  labs(title = 'True Second-Order ALEs (excluding Main Effects)',
       x = 'Predictor 1 Value',
       y = 'Predictor 2 Value',
       fill = 'ALE')

# Estimate (with Main Effects)
limits <- c(min(ale2$truth_main, ale2$est_main), max(ale2$truth_main, ale2$est_main))
ale2_est_main <- ale2 |>
  ggplot(aes(xmin = x1 - w1, xmax = x1 + w2, ymin = x2 - h1, ymax = x2 + h2)) +
  geom_rect(aes(fill = est_main)) +
  theme_bw() +
  coord_fixed() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_gradient2(low = 'blue', high = 'red', limits = limits) +
  facet_grid(var2 ~ var1) +
  theme_bw() +
  labs(title = 'Estimated Second-Order ALEs (including Main Effects)',
       x = 'Predictor 1 Value',
       y = 'Predictor 2 Value',
       fill = 'ALE')

# Truth (with Main Effects)
ale2_true_main <- ale2 |>
  ggplot(aes(xmin = x1 - w1, xmax = x1 + w2, ymin = x2 - h1, ymax = x2 + h2)) +
  geom_rect(aes(fill = truth_main)) +
  theme_bw() +
  coord_fixed() +
  scale_x_continuous(expand = c(0,0), breaks = seq(0,1,0.2)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_gradient2(low = 'blue', high = 'red', limits = limits) +
  facet_grid(var2 ~ var1) +
  theme_bw() +
  labs(title = 'True Second-Order ALEs (including Main Effects)',
       x = 'Predictor 1 Value',
       y = 'Predictor 2 Value',
       fill = 'ALE')

library(patchwork)
ale2_est + ale2_true
ale2_est_main + ale2_true_main
