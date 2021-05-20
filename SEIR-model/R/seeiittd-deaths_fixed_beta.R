library(cmdstanr)
library(dplyr)
library(lubridate)
library(LaplacesDemon)

# Fixed Beta SEEIITTD
model <- cmdstan_model(here::here("SEIR-model", "stan", "deaths_fixed_beta_seeiittd.stan"))

# Real data
br <- readRDS(here::here("SEIR-model/", "data", "brazil_nation.rds"))

initial_time <- 0
total_time <- br %>% nrow
n_disease_states = 8 # seeiittd
population <- br %>% pull(estimated_population_2019) %>% max

#Deaths are reported weekly, reporting periods are every seven days
deaths_starts <- seq(1, (total_time / 7 - 1) * 7, by = 7)
deaths_stops <- seq(7, total_time, by = 7)
deaths_length <- min(c(length(deaths_starts), length(deaths_starts)))
deaths <- br %>%
  group_by(week = cut(date, "week")) %>%
  summarise(deaths = sum(new_deaths)) %>%
  pull(deaths) %>%
  ceiling %>%
  as.integer %>% 
  head(-1)

stan_data <- list(
  initial_time = initial_time,
  T = total_time,
  times = seq_len(total_time),
  n_disease_states = n_disease_states,
  population = population,
  deaths_length = deaths_length,
  deaths_starts = deaths_starts,
  deaths_stops = deaths_stops,
  deaths = deaths
)

fit <- model$sample(data = stan_data,
                        seed = 321,
                        init = function() list(initial_state_raw = c(runif(1, min = 0.99999, max = 1.0), runif(1, min = 0.0, max = 1.0)),
                                               beta = exp(runif(1, min = -2, max = 0.5)),
                                               dL = runif(1, min = 3.5, max = 4.0),
                                               dI = runif(1, min = 2.2, max = 2.6),
                                               dT = runif(1, min = 11.0, max = 13.0),
                                               omega = invlogit(runif(1, min = -5, max = -3))),
                    chains = 4,
                    parallel_chains = 4,
                    output_dir = here::here("SEIR-model", "results", "deaths_fixed_beta", "seeiittd"))
