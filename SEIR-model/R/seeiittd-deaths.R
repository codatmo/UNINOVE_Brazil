library(cmdstanr)
library(dplyr)
library(rstan)

# Simulated Data from Liverpool's CoDatMo
# src/scripts/generateSimulatedData.R
# generatedData <- generateSimulatedData(
#   seed = 0,
#   maxTime = 100,
#   population = 1000000,
#   n_beta_pieces = 20,
#   n_rho_calls_111_pieces = 10,
#   calls_111_start = 30
# )
# stan_data_sim <- generatedData$stan_data

stan_data_sim <- readRDS(here::here("SEIR-model", "data", "stan_data_sim.rds"))
stan_data_sim$compute_likelihood <- 1
stan_data_sim$n_rho_calls_111_pieces <- NULL
stan_data_sim$rho_calls_111_left_t <- NULL
stan_data_sim$rho_calls_111_right_t <- NULL
stan_data_sim$calls_111_length <- NULL
stan_data_sim$calls_111_start <- NULL
stan_data_sim$calls_111 <- NULL

model <- cmdstan_model(here::here("SEIR-model", "stan", "deaths.stan"))
fit_sim <- model$sample(data = stan_data_sim,
                        seed = 123,
                        chains = 4)

summary_model_sim <- fit_sim$summary()

summary_model_sim %>%
  filter(rhat >= 1.01)



fit <- model$sample(data = stan_data,
                        seed = 123,
                        chains = 4,
                        output_dir = here::here("SEIR-model", "results", "deaths"))

# If necessary you can load with
files <- list.files(here::here("SEIR-model", "results", "deaths"), full.names = TRUE)
output <- read_cmdstan_csv(files)

# Plotting Stuff
r_stan_sir <- rstan::read_stan_csv(files)