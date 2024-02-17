
library(tidyverse)
library(patchwork)

results <- read_rds(here::here('Results','16FEB2024.rds'))

params <- read_csv(here::here('Params','params-new.csv'))

sim_stats <- lapply(results, '[[', 'stats') |>
  bind_rows(.id = 'key') |>
  mutate(key = as.numeric(key)) |>
  left_join(params, by = join_by('key' == 'key'))

# Bias
sim_stats |>
  filter(!(param %in% c('G','nu'))) |>
  ggplot(aes(x = factor(m), y = bias_1, ymin = bias_2, ymax = bias_3, 
             color = soft, lty = sparse, shape = sparse)) +
  geom_pointrange(size = 0.3, position = position_dodge2(0.5)) + 
  facet_wrap(~param, scales = 'free', nrow = 2) +
  theme_bw() +
  theme(legend.position = 'bottom') +
  labs(x = 'Number of Trees',
       y = 'Bias',
       color = 'Soft Trees?',
       lty = 'Sparse Trees?',
       shape = 'Sparse Trees?')

# Coverage
sim_stats |>
  filter(!(param %in% c('G','nu'))) |>
  ggplot(aes(x = factor(m), y = coverage_1, color = soft, shape = sparse)) +
  geom_hline(yintercept = 0.95, lty = 2) +
  geom_point() +
  facet_wrap(~param, scales = 'free', nrow = 2) +
  theme_bw() +
  theme(legend.position = 'bottom') +
  labs(x = 'Number of Trees',
       y = '95% Credible Interval Coverage',
       color = 'Soft Trees?',
       lty = 'Sparse Trees?',
       shape = 'Sparse Trees?')

# BART predictions
sim_stats |>
  filter(param == 'G') |>
  pivot_longer(cols = bias_1:rmse_3,
               names_to = c(".value", "number"),
               names_pattern = "([a-z]+)_(\\d+)") |>
  pivot_longer(cols = bias:rmse, names_to = 'type') |>
  pivot_wider(id_cols = c(key, type, m, sparse, soft), names_from = number, values_from = value) |>
  rename(est = `1`, lower = `2`, upper = `3`) |>
  ggplot(aes(x = factor(m), y = est, ymin = lower, ymax = upper, color = soft, lty = sparse, shape = sparse)) +
  geom_pointrange(size = 0.3, position = position_dodge2(0.5)) +
  facet_wrap(~type, scales = 'free') +
  theme_bw() +
  theme(legend.position = 'bottom',
        text = element_text(size = 13)) +
  labs(x = 'Number of Trees',
       y = '95% Credible Interval Coverage',
       color = 'Soft Trees?',
       lty = 'Sparse Trees?',
       shape = 'Sparse Trees?')

# First order ALE
results[[12]]$ale1 |>
  filter(!(x %in% c(min(x), max(x))), .by = var) |>
  ggplot(aes(x = x, ymin = lcl, ymax = ucl)) +
  geom_hline(yintercept = 0, lty = 2, color = 'gray') +
  geom_ribbon(alpha = 0.2) +
  geom_line(aes(y = est, col = 'Estimate')) +
  geom_line(aes(y = truth, col = 'Truth'), lty = 2) +
  facet_wrap(~factor(var, levels = unique(ale1$var)), ncol = 5, scales = 'free_y') +
  theme_bw() +
  theme(legend.position = 'bottom',
        text = element_text(size = 14)) +
  scale_color_manual(values = c("Truth" = "red", "Estimate" = "black")) +
  scale_x_continuous(breaks = c(0.25, 0.5, 0.75)) +
  labs(#title = 'ALE Plots of First-Order Main Effects',
       x = 'Predictor Value',
       y = 'ALE Main Effect (95% Credible Interval)',
       color = '') +
  guides(color = guide_legend(override.aes = list(lty = c(1, 2))))

# Second order ALE
lapply(results[12], \(x) {
  ale2 <- x$ale2 |>
    filter(!(x1 %in% c(min(x1), max(x1)) | x2 %in% c(min(x2), max(x2))), .by = c(var1, var2))
  limits <- c(min(ale2$truth, ale2$est), max(ale2$truth, ale2$est))
  
  est <- ale2 |>
    ggplot(aes(xmin = x1 - w1, xmax = x1 + w2, ymin = x2 - h1, ymax = x2 + h2)) +
    geom_rect(aes(fill = est), show.legend = FALSE) +
    theme_bw() +
    coord_fixed() +
    scale_x_continuous(expand = c(0,0), breaks = c(0.25, 0.5, 0.75)) +
    scale_y_continuous(expand = c(0,0), breaks = c(0.25, 0.5, 0.75)) +
    scale_fill_gradient2(low = 'blue', high = 'red', limits = limits) +
    facet_grid(var2 ~ var1) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 270),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          text = element_text(size = 13)) +
    labs(# title = 'Estimated Second-Order ALEs',
         #subtitle = 'Excluding Main Effects',
      title = '(a)',
         x = 'Predictor 1 Value',
         y = 'Predictor 2 Value',
         fill = 'ALE')
  
  true <- ale2 |>
    ggplot(aes(xmin = x1 - w1, xmax = x1 + w2, ymin = x2 - h1, ymax = x2 + h2)) +
    geom_rect(aes(fill = truth)) +
    theme_bw() +
    coord_fixed() +
    scale_x_continuous(expand = c(0,0), breaks = c(0.25, 0.5, 0.75)) +
    scale_y_continuous(expand = c(0,0), breaks = c(0.25, 0.5, 0.75)) +
    scale_fill_gradient2(low = 'blue', high = 'red', limits = limits) +
    facet_grid(var2 ~ var1) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 270),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          text = element_text(size = 13)) +
    labs(#title = 'True Second-Order ALEs',
         #subtitle = 'Excluding Main Effects',
          title = '(b)',
         x = 'Predictor 1 Value1',
         y = 'Predictor 2 Value',
         fill = 'ALE')
  
  est + true
})


# First + Second order ALE
lapply(results[12], \(x) {
  ale3 <- x$ale3 |>
    filter(!(x1 %in% c(min(x1), max(x1)) | x2 %in% c(min(x2), max(x2))), .by = c(var1, var2))
  limits <- c(min(ale3$truth, ale3$est), max(ale3$truth, ale3$est))
  
  est <- ale3 |>
    ggplot(aes(xmin = x1 - w1, xmax = x1 + w2, ymin = x2 - h1, ymax = x2 + h2)) +
    geom_rect(aes(fill = est), show.legend = FALSE) +
    theme_bw() +
    coord_fixed() +
    scale_x_continuous(expand = c(0,0), breaks = c(0.25, 0.5, 0.75)) +
    scale_y_continuous(expand = c(0,0), breaks = c(0.25, 0.5, 0.75)) +
    scale_fill_gradient2(low = 'blue', high = 'red', limits = limits) +
    facet_grid(var2 ~ var1) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 270),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          text = element_text(size = 13)) +
    labs(#title = 'True Second-Order ALEs',
      #subtitle = 'Excluding Main Effects',
      title = '(a)',
      x = 'Predictor 1 Value1',
      y = 'Predictor 2 Value',
      fill = 'ALE')
  
  true <- ale3 |>
    ggplot(aes(xmin = x1 - w1, xmax = x1 + w2, ymin = x2 - h1, ymax = x2 + h2)) +
    geom_rect(aes(fill = truth)) +
    theme_bw() +
    coord_fixed() +
    scale_x_continuous(expand = c(0,0), breaks = c(0.25, 0.5, 0.75)) +
    scale_y_continuous(expand = c(0,0), breaks = c(0.25, 0.5, 0.75)) +
    scale_fill_gradient2(low = 'blue', high = 'red', limits = limits) +
    facet_grid(var2 ~ var1) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 270),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          text = element_text(size = 13)) +
    labs(#title = 'True Second-Order ALEs',
      #subtitle = 'Excluding Main Effects',
      title = '(b)',
      x = 'Predictor 1 Value1',
      y = 'Predictor 2 Value',
      fill = 'ALE')
  
  est + true
})
