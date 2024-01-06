# Program Name: spanbbart-sim-analysis.R
# Author:       Jacob Englert
# Date:         2 January 2024
# Purpose:      Analyze results of the spanbbart simulation.


# Load Packages -----------------------------------------------------------
library(tidyverse)
library(patchwork)


# Load Results ------------------------------------------------------------
sim_date <- '05JAN2024'
results <- read_rds(here::here('Results', paste0(sim_date, '.rds')))


# Summarize Simulation Statistics -----------------------------------------

# Bias, coverage, and RMSE
lower <- function(x) quantile(x, 0.025, na.rm = TRUE)
upper <- function(x) quantile(x, 0.975, na.rm = TRUE)
stats <- lapply(results, '[[', 'stats') |>
  bind_rows() |>
  group_by(param) |>
  summarise(across(c(bias, coverage, rmse), list(mean, lower, upper)))


# Summarize First Order ALE -----------------------------------------------

# First-order ALE
ale1 <- lapply(results, '[[', 'ale1') |>
  bind_rows() |>
  summarise(across(c(est, lcl, ucl, truth), mean), .by = c(x, var))

ale1 |>
  filter(!(x %in% c(min(x), max(x))), .by = var) |>
  ggplot(aes(x = x, ymin = lcl, ymax = ucl)) +
  geom_hline(yintercept = 0, lty = 2, color = 'gray') +
  geom_line(aes(y = est, col = 'Estimate')) +
  geom_line(aes(y = truth, col = 'Truth'), lty = 2) +
  geom_ribbon(alpha = 0.4) +
  facet_wrap(~factor(var, levels = unique(ale1$var)), ncol = 5, scales = 'free_y') +
  theme_bw() +
  theme(legend.position = 'bottom') +
  scale_color_manual(values = c("Truth" = "red", "Estimate" = "black")) +
  labs(title = 'ALE Plots of First-Order Main Effects',
       x = 'Predictor Value',
       y = 'ALE',
       color = 'Type') +
  guides(color = guide_legend(override.aes = list(lty = c(1, 2))))
ggsave(here::here('Figures', sim_date,'ale-first-order.png'), width = 6.5, height = 5)


# Summarize Second Order ALE ----------------------------------------------

# Second-order ALE
ale2 <- lapply(results, '[[', 'ale2') |>
  bind_rows() |>
  summarise(across(c(est, lcl, ucl, truth, est_main, truth_main), mean), 
            .by = c(x1, x2, w1, w2, h1, h2, var1, var2)) |>
  arrange(var1, var2, x1, x2) |>
  filter(!(x1 %in% c(min(x1), max(x1)) | x2 %in% c(min(x2), max(x2))), .by = c(var1, var2))


# Estimate (without Main Effects)
limits <- c(min(ale2$truth, ale2$est), max(ale2$truth, ale2$est))
ale2_est <- ale2 |>
  ggplot(aes(xmin = x1 - w1, xmax = x1 + w2, ymin = x2 - h1, ymax = x2 + h2)) +
  geom_rect(aes(fill = est), show.legend = FALSE) +
  theme_bw() +
  coord_fixed() +
  scale_x_continuous(expand = c(0,0), breaks = seq(0,1,0.25)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_gradient2(low = 'blue', high = 'red', limits = limits) +
  facet_grid(var2 ~ var1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 270),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  labs(title = 'Estimated Second-Order ALEs',
       subtitle = 'Excluding Main Effects',
       x = 'Predictor 1 Value',
       y = 'Predictor 2 Value',
       fill = 'ALE')

# Truth (without Main Effects)
ale2_true <- ale2 |>
  ggplot(aes(xmin = x1 - w1, xmax = x1 + w2, ymin = x2 - h1, ymax = x2 + h2)) +
  geom_rect(aes(fill = truth)) +
  theme_bw() +
  coord_fixed() +
  scale_x_continuous(expand = c(0,0), breaks = seq(0,1,0.25)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_gradient2(low = 'blue', high = 'red', limits = limits) +
  facet_grid(var2 ~ var1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 270),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  labs(title = 'True Second-Order ALEs',
       subtitle = 'Excluding Main Effects',
       x = 'Predictor 1 Value',
       y = 'Predictor 2 Value',
       fill = 'ALE')

ale2_est + ale2_true + plot_layout(guides = 'collect')
ggsave(here::here('Figures', sim_date,'ale-second-order.png'), width = 6.5, height = 4)


# Estimate (with Main Effects)
limits <- c(min(ale2$truth_main, ale2$est_main), max(ale2$truth_main, ale2$est_main))
ale2_est_main <- ale2 |>
  ggplot(aes(xmin = x1 - w1, xmax = x1 + w2, ymin = x2 - h1, ymax = x2 + h2)) +
  geom_rect(aes(fill = est_main)) +
  theme_bw() +
  coord_fixed() +
  scale_x_continuous(expand = c(0,0), breaks = seq(0,1,0.25)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_gradient2(low = 'blue', high = 'red', limits = limits) +
  facet_grid(var2 ~ var1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 270),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  labs(title = 'Estimated Second-Order ALEs',
       subtitle = 'Including Main Effects',
       x = 'Predictor 1 Value',
       y = 'Predictor 2 Value',
       fill = 'ALE')

# Truth (with Main Effects)
ale2_true_main <- ale2 |>
  ggplot(aes(xmin = x1 - w1, xmax = x1 + w2, ymin = x2 - h1, ymax = x2 + h2)) +
  geom_rect(aes(fill = truth_main)) +
  theme_bw() +
  coord_fixed() +
  scale_x_continuous(expand = c(0,0), breaks = seq(0,1,0.25)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_gradient2(low = 'blue', high = 'red', limits = limits) +
  facet_grid(var2 ~ var1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 270),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  labs(title = 'True Second-Order ALEs',
       subtitle = 'Including Main Effects',
       x = 'Predictor 1 Value',
       y = 'Predictor 2 Value',
       fill = 'ALE')

# Plot posterior mean and truth side-by-side
ale2_est_main + ale2_true_main + plot_layout(guides = 'collect')
ggsave(here::here('Figures', sim_date,'ale-second-order-w-main.png'), width = 6.5, height = 4)
