library(cmdstanr)
library(dplyr)
library(lubridate)
#library(LaplacesDemon)

# RK4 Custom Model(https://github.com/spinkney/helpful_stan_functions/blob/main/functions/ode/odeint_rk4.stan)
# based on Funko_Unko's contribution: https://discourse.mc-stan.org/t/codatmo-liverpool-uninove-models-slow-ode-implementation-and-trapezoidal-solver/22500/42
model <- cmdstan_model(here::here("SEIR-model", "stan", "deaths_penalized_rk4.stan"))

# The interface is mostly the same as before (rk4 model). He did add some things:
#   
# int beta_linear_discontinuous;
# int beta_linear_continuous;
# int beta_constant_discontinuous;
# real beta_regularization;
# int num_sub_steps;
# 
# Setting to one either of `beta_linear_discontinuous`, `beta_linear_continuous` or `beta_constant_discontinuous` selects the corresponding ansatz for beta
# `num_sub_steps` lets you perform intermediate `rk4` steps. Set to one for the original behavior
# For beta_regularization: there are three behaviors:
#   
# 1. beta_regularization = 0: no regularization
# 2. beta_regularization > 0: do something like beta_left[2:] ~ normal(beta_left[:n_beta_pieces-1], dummy_beta_reg);
# 3. beta_regularization < 0: in addition treat beta_regularization as a log-normally distributed parameter.
# 
# However, everything is in serious need of refactoring.
# 
# For fitting, specifying metric=‘dense_e’ appears to speed things up and I think inits should be samples from the prior with constant-in-time betas.

# Real data
br <- readRDS(here::here("SEIR-model/", "data", "brazil_nation.rds"))

initial_time <- 0
total_time <- br %>% nrow
beta_pieces <- 7 # number of days to make beta_pieces (weekly)
n_beta_pieces <- ceiling(total_time / beta_pieces)
beta_left_t <- seq(0, total_time - 1, by = beta_pieces)
beta_right_t <- c(beta_left_t[2:n_beta_pieces], total_time + 1)
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

# I have no idea what real data is...
real_data_length <- length(beta_left_t) + length(beta_right_t) + 1
real_data <- c(beta_left_t, beta_right_t, population)
integer_data_length <- 5
integer_data <- c(total_time, length(beta_left_t), length(beta_right_t), length(beta_left_t), n_disease_states)

stan_data <- list(
  initial_time = initial_time,
  n_beta_pieces = n_beta_pieces,
  beta_left_t = beta_left_t,
  beta_right_t = beta_right_t,
  T = total_time,
  times = seq_len(total_time),
  n_disease_states = 8,
  population = population,
  deaths_length = deaths_length,
  deaths_starts = deaths_starts,
  deaths_stops = deaths_stops,
  deaths = deaths,
  real_data_length = real_data_length,
  real_data = real_data,
  integer_data_length = integer_data_length,
  integer_data = integer_data,
  compute_likelihood = 1,
  beta_linear_discontinuous = 1,
  beta_linear_continuous = 0,
  beta_constant_discontinuous = 0,
  beta_regularization = 0.2, # <<< Most Important!
  num_sub_steps = 1 # how many betas in weeks will be used?
)

fit <- model$sample(data = stan_data,
                    seed = 321,
                    metric="dense_e",
                    # init = function() list(initial_state_raw = c(runif(1, min = 0.99999, max = 1.0), runif(1, min = 0.0, max = 1.0)),
                    #                        beta_left = exp(runif(n_beta_pieces, min = -2, max = 0.5)),
                    #                        beta_right = exp(runif(n_beta_pieces, min = -2, max = 0.5)),
                    #                        dL = runif(1, min = 3.5, max = 4.0),
                    #                        dI = runif(1, min = 2.2, max = 2.6),
                    #                        dT = runif(1, min = 11.0, max = 13.0),
                    #                        omega = invlogit(runif(1, min = -5, max = -3))),
                    chains = 4,
                    parallel_chains = 4,
                    output_dir = here::here("SEIR-model", "results", "deaths_penalized_rk4"))

# If necessary you can load with
# files <- list.files(here::here("SEIR-model", "results", "deaths_rk4"), full.names = TRUE)
# output <- read_cmdstan_csv(files)
