library(rstan)
library(cmdstanr)
library(dplyr)
library(tibble)

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

stan_data_sim_deaths <- readRDS(here::here("SEIR-model", "data", "stan_data_sim_deaths.rds"))
stan_data_sim_twitter <- readRDS(here::here("SEIR-model", "data", "stan_data_sim_twitter.rds"))

# Trapeizodal Rule Model
model_deaths <- cmdstan_model(here::here("SEIR-model", "stan", "deaths_trapezoidal.stan"))
model_twitter <- cmdstan_model(here::here("SEIR-model", "stan", "deaths_twitter_trapezoidal.stan"))

# Run Models
fit_sim_deaths <- model_deaths$sample(data = stan_data_sim_deaths,
                                        seed = 123,
                                        parallel_chains = 4)

fit_sim_twitter <- model_twitter$sample(data = stan_data_sim_twitter,
                                        seed = 123,
                                        parallel_chains = 4)

# Compare MAE Predicted vs Real Deaths
real_deaths <- stan_data_sim_deaths$deaths %>%
  enframe(name = "day", value = "real_deaths")

fit_sim_deaths$summary("pred_deaths") %>% 
  bind_cols(real_deaths) %>% 
  mutate(
    MAE_median = abs(median - real_deaths),
    MAE_mean = abs(mean - real_deaths)) %>% 
  summarise(
    MAE_median = mean(MAE_median),
    MAE_mean = mean(MAE_mean))

fit_sim_twitter$summary("pred_deaths") %>% 
  bind_cols(real_deaths) %>% 
  mutate(
    MAE_median = abs(median - real_deaths),
    MAE_mean = abs(mean - real_deaths)) %>% 
  summarise(
    MAE_median = mean(MAE_median),
    MAE_mean = mean(MAE_mean))
