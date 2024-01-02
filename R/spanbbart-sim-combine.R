#!/usr/bin/env Rscript

library(tidyverse)

today <- toupper(format(Sys.Date(), "%d%b%Y"))

results <- list.files(here::here('Results','temp'), full.names = TRUE) |>
  lapply(read_rds)

write_rds(results, here::here('Results', paste0(today, '.rds')))