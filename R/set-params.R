# Program Name: set-params.R
# Author: Jacob Englert
# Date: 02JAN2024
# Description: Generate parameter tables for spanbbart simulations.

# Load Packages -----------------------------------------------------------
library(tidyverse)


# Specify Simulation Parameters -------------------------------------------
nT <- 300
m <- 200
k <- 2
power <- 2
base <- 0.95
K <- 100
seed <- 1:500

params <- crossing(nT, m, k, power, base, K, seed) |>
  mutate(num_iter = 5000, num_burn = 2500, num_thin = 5,
         ID = row_number()) |>
  select(ID, everything())

write_csv(params, here::here('Params','params.csv'))
