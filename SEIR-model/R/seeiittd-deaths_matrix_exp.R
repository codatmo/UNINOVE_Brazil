library(cmdstanr)
library(dplyr)
library(lubridate)


# based on Funko_Unko's contribution: https://discourse.mc-stan.org/t/codatmo-liverpool-uninove-models-slow-ode-implementation-and-trapezoidal-solver/22500/45
model <- cmdstan_model(here::here("SEIR-model", "stan", "deaths_matrix_exp.stan"))

# Real data
br <- readRDS(here::here("SEIR-model/", "data", "brazil_nation.rds"))

no_days <- br %>% nrow
population <- br %>% pull(estimated_population_2019) %>% max
deaths <- br %>% pull(new_deaths) %>% cumsum

stan_data <- list(
  no_days = no_days,
  population = population,
  deaths = deaths
)

fit <- model$sample(data = stan_data,
                    seed = 321,
                    parallel_chains = 4,
                    output_dir = here::here("SEIR-model", "results", "deaths_matrix_exp"))

# If necessary you can load with
# files <- list.files(here::here("SEIR-model", "results", "deaths_matrix_exp"), full.names = TRUE)
# output <- read_cmdstan_csv(files)
