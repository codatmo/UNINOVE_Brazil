library(rstan)
library(cmdstanr)
library(dplyr)
library(tibble)

source(here::here("SEIR-model", "R", "generateSimulatedData.R"))

stan_data <- generateSimulatedData(seed = 123, n_beta_pieces = 20, maxTime = 100)$stan_data
ground_truth <- generateSimulatedData(seed = 123, n_beta_pieces = 20, maxTime = 100)$ground_truth
stan_data$compute_likelihood <- 1
stan_data$compute_twitter <- 1

# Trapeizodal Rule Model
model_deaths <- cmdstan_model(here::here("SEIR-model", "stan", "deaths_trapezoidal.stan"))
model_twitter <- cmdstan_model(here::here("SEIR-model", "stan", "deaths_twitter_trapezoidal_informative.stan"))

# Run Models
fit_sim_deaths <- model_deaths$sample(data = stan_data,
                                        seed = 123,
                                        parallel_chains = 4)

fit_sim_twitter <- model_twitter$sample(data = stan_data,
                                        seed = 123,
                                        parallel_chains = 4)

# Compare MAE Predicted vs Real Deaths
real_deaths <- stan_data$deaths %>%
  enframe(name = "day", value = "real_deaths")

fit_sim_deaths$summary("pred_deaths") %>%
  bind_cols(real_deaths) %>%
  mutate(
    MAE_median = abs(median - real_deaths),
    MAE_mean = abs(mean - real_deaths)) %>%
  summarise(
    MAE_median = mean(MAE_median),
    MAE_mean = mean(MAE_mean)
  )


fit_sim_twitter$summary("pred_deaths") %>%
  bind_cols(real_deaths) %>%
  mutate(
    MAE_median = abs(median - real_deaths),
    MAE_mean = abs(mean - real_deaths)) %>%
  summarise(
    MAE_median = mean(MAE_median),
    MAE_mean = mean(MAE_mean)
  )

# Compare dI
fit_sim_deaths$summary("dI")
fit_sim_twitter$summary("dI")

# Compare dT
fit_sim_deaths$summary("dT")
fit_sim_twitter$summary("dT")
