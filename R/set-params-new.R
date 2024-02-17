# Program Name: set-params.R
# Author: Jacob Englert
# Date: 02JAN2024
# Description: Generate parameter tables for spanbbart simulations.

# Load Packages -----------------------------------------------------------
library(tidyverse)


# Specify Simulation Parameters -------------------------------------------
nT <- 300
m <- c(10, 25, 50, 100)
k <- 2
power <- 2
base <- 0.95
K <- 40
soft <- c(TRUE, FALSE)
sparse <- c(TRUE, FALSE)
num_iter <- 10000
num_burn <- 5000
num_thin <- 10

params <- crossing(nT, m, k, power, base, K, soft, sparse,
                   num_iter, num_burn, num_thin) |>
  mutate(key = row_number()) |>
  select(key, everything())

write_csv(params, here::here('Params','params-new.csv'))
