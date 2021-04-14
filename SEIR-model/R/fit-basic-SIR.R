library(cmdstanr)
library(dplyr)

br <- readRDS(here::here("SEIR-model/", "data", "brazil_2021.rds"))

# total count
N <- br %>% pull(estimated_population_2019) %>% max

# time series of cases
cases <- br %>% pull(new_confirmed)

# times
n_days <- br %>% pull(date) %>% unique %>% length
ts <- n_days %>% seq_len
t0 <- 0

#initial conditions
i0 <- br %>% pull(new_confirmed) %>% min # Infected
s0 <- N - i0 # Susceptible
r0 <- 0 # Resolved
y0 = c(S = s0, I = i0, R = r0)

compute_likelihood = 1

# data for Stan
data_sir <- list(n_days = n_days, i0 = i0, y0 = y0, s0 = s0, r0 = r0, t0 = t0, ts = ts, 
                 N = N, cases = cases, compute_likelihood = compute_likelihood)

# Compile model
model <- cmdstan_model(here::here("SEIR-model", "stan", "basic-SIR.stan"))

# Fit model
fit_sir <- model$sample(data = data_sir)

# Summary Stats
fit_sir$summary()
fit_sir$summary("R0")

# Save the model as CSV files
fit_sir$save_output_files(here::here("SEIR-model", "results", "basic-SIR"))

# If necessary you can load with
files <- list.files(here::here("SEIR-model", "results", "basic-SIR"), full.names = TRUE)
output <- read_cmdstan_csv(files)
