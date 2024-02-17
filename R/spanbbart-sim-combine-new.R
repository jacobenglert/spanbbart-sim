#!/usr/bin/env Rscript

library(tidyverse)


# All result files
all_files <- list.files(here::here('Results','temp'))

# Allocate storage for summarized results
results <- list()

# Keys included in results
keys <- sort(unique(as.numeric(substr(all_files, 1, 3))))

# Iterate through each key
for (key in keys) {
  
  # Read key files
  key_files <- all_files[which(as.numeric(substr(all_files, 1, 3)) == key)]
  key_results <- lapply(key_files, \(f) read_rds(here::here('Results','temp', f)))
  
  # Summarize bias, coverage, and RMSE
  lower <- function(x) quantile(x, 0.025, na.rm = TRUE)
  upper <- function(x) quantile(x, 0.975, na.rm = TRUE)
  stats <- lapply(key_results, '[[', 'stats') |>
    bind_rows() |>
    group_by(param) |>
    summarise(across(c(est, bias, coverage, rmse), list(mean, lower, upper)))
  
  # Summarize first order ALE
  ale1 <- lapply(key_results, '[[', 'ale1') |>
    bind_rows() |>
    summarise(across(c(est, lcl, ucl, truth), mean), .by = c(x, var))

  # Summarize second order ALE
  ale2 <- lapply(key_results, '[[', 'ale2') |>
    bind_rows() |>
    summarise(across(c(est, lcl, ucl, truth), mean), 
              .by = c(x1, x2, w1, w2, h1, h2, var1, var2)) |>
    arrange(var1, var2, x1, x2)
  
  # Summarize first + second order ALE
  ale3 <- lapply(key_results, '[[', 'ale3') |>
    bind_rows() |>
    summarise(across(c(est, lcl, ucl, truth), mean), 
              .by = c(x1, x2, w1, w2, h1, h2, var1, var2)) |>
    arrange(var1, var2, x1, x2)
  
  results[[key]] <- list(key = key, stats = stats, ale1 = ale1, ale2 = ale2, ale3 = ale3)
  
}

# Write combined results file
today <- toupper(format(Sys.Date(), "%d%b%Y"))
write_rds(results, here::here('Results', paste0(today, '.rds')))