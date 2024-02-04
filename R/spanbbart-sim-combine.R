#!/usr/bin/env Rscript

library(tidyverse)

# Today's date
today <- toupper(format(Sys.Date(), "%d%b%Y"))

# Result files
results_files <- list.files(here::here('Results','temp'), full.names = TRUE)

# Read result files
results <- lapply(results_files, read_rds)

# Delete result files
lapply(results_files, file.remove)
lapply(list.files(here::here('Results','temp'), full.names = TRUE), file.remove)

# Write combined results file
write_rds(results, here::here('Results', paste0(today, '.rds')))