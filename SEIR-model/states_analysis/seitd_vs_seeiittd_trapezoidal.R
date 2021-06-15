library(cmdstanr)
library(dplyr)
library(lubridate)
library(readr)
library(tibble)
library(tidyr)
library(rstan)

seitd <- list.files(here::here("SEIR-model", "results", "deaths_trapezoidal", "seitd"), full.names = TRUE) %>%
    as_cmdstan_fit()

dt_setid <- seitd$summary(variables = "dT")
seitd <- seitd$summary(variables = c("S", "E", "I", "T_", "D"))

seeiittd <- list.files(here::here("SEIR-model", "results", "deaths_trapezoidal", "seeiittd"), full.names = TRUE) %>%
    as_cmdstan_fit()

dt_seeiittd <- seeiittd$summary(variables = "dT")
seeiittd <- seeiittd$summary(variables = c("S",
                                           "E1", "E2",
                                           "I1", "I2",
                                           "T1", "T2",
                                           "D"))

seitd <- seitd %>%
    extract(variable, c("state", "day"),
            "state_(\\w)\\[([\\d{1}]+)\\]",
            convert = TRUE) %>%
            select(state, day, median) %>%
            pivot_wider(names_from = state,
                        values_from = median) %>%
            mutate(day = seq(ymd("2020-02-25"), ymd("2021-05-04"), by = "1 day"))

seeiittd <- seeiittd %>%
    extract(variable, c("state", "day"),
            "state_(\\w[\\d]?)\\[([\\d{1}]+)\\]",
            convert = TRUE) %>%
            select(state, day, median) %>%
            pivot_wider(names_from = state,
                        values_from = median) %>%
            mutate(day = seq(ymd("2020-02-25"), ymd("2021-05-04"), by = "1 day"))

# Real Deaths
br <- readRDS(here::here("SEIR-model/", "data", "brazil_nation_2020.rds"))
real_deaths <- br %>%
    pull(new_deaths) %>%
    cumsum %>%
    enframe(name = "day", value = "real_deaths") %>%
    mutate(day = seq(ymd('2020-02-25'),ymd('2021-05-04'), by = '1 day'))

real_deaths %>%
    left_join(seitd, "day") %>%
    left_join(seeiittd, "day",
    suffix = c("_seitd", "_seeiitd")) %>%
    write_csv(here::here("SEIR-model", "states_analysis", "trapezoidal_seitd_vs_seeiittd.csv"))

# Plot Stuff
rstan_seitd <- rstan::read_stan_csv(list.files(here::here("SEIR-model",
                                                          "results",
                                                          "deaths_trapezoidal",
                                                          "seitd"),
                                               full.names = TRUE))
rstan_seeiittd <- rstan::read_stan_csv(list.files(here::here("SEIR-model",
                                                             "results",
                                                             "deaths_trapezoidal",
                                                             "seeiittd"),
                                                  full.names = TRUE))
stan_dens(rstan_seitd, pars = "dT", separate_chains = FALSE)
stan_dens(rstan_seeiittd, pars = "dT", separate_chains = FALSE)
