library(cmdstanr)
library(dplyr)
library(lubridate)
library(readr)
library(tibble)
library(tidyr)
library(rstan)
library(loo)

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

seitd_wide <- seitd %>%
    tidyr::extract(variable, c("state", "day"),
            "(.*)\\[([\\d{1}]+)\\]",
            convert = TRUE) %>%
    select(state, day, median) %>%
    pivot_wider(names_from = state,
                        values_from = median) %>%
    mutate(day = seq(ymd("2020-02-25"), ymd("2021-01-01"), by = "1 day")) %>%
    rename(T = T_) %>%
    head(-1)

seitd_long <- seitd_wide %>%
    pivot_longer(-day, names_to = "state", values_to = "median")

seeiittd_wide <- seeiittd %>%
    tidyr::extract(variable, c("state", "day"),
            "(\\w[\\d]?)\\[([\\d{1}]+)\\]",
            convert = TRUE) %>%
    select(state, day, median) %>%
    pivot_wider(names_from = state,
                values_from = median) %>%
    mutate(day = seq(ymd("2020-02-25"), ymd("2021-01-01"), by = "1 day")) %>%
    head(-1)  %>%
    transmute(
        day, S,
        E = E1 + E2,
        I = I1 + I2,
        T = T1 + T2,
        D
    )

seeiittd_long <- seeiittd_wide %>%
    pivot_longer(-day, names_to = "state", values_to = "median")

# Real Deaths
br <- readRDS(here::here("SEIR-model", "data", "brazil_nation_2020.rds"))
real_deaths <- br %>%
    pull(new_deaths) %>%
    cumsum %>%
    enframe(name = "day", value = "real_deaths") %>%
    mutate(day = seq(ymd('2020-02-25'),ymd('2020-12-31'), by = '1 day'))

real_deaths %>%
    left_join(seitd, "day") %>%
    left_join(seeiittd, "day",
    suffix = c("_seitd", "_seeiitd")) %>%
    write_csv(here::here("SEIR-model", "states_analysis", "trapezoidal_seitd_vs_seeiittd.csv"))

# Plot Stuff

seitd_long %>% mutate(group = "seitd") %>%
    bind_rows(seeiittd_long  %>% mutate(group = "seeiittd"))  %>%
    filter(state != "S") %>%
    ggplot(aes(x = day, y = median, color = state, linetype = group)) +
    geom_line() +
    scale_color_brewer(palette = "Set1") +
    scale_y_continuous(labels = scales::label_number_si()) +
    theme(legend.position = "bottom")

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
stan_dens(rstan_seitd,
          pars = "dT",
          separate_chains = FALSE) +
          ggtitle("SEITD")
stan_dens(rstan_seeiittd,
          pars = "dT",
          separate_chains = FALSE) +
          ggtitle("SEEIITTD")

# LOO-CV
# Extract pointwise log-likelihood
# using merge_chains=FALSE returns an array, which is easier to
# use with relative_eff()
log_lik_seitd <- extract_log_lik(rstan_seitd, merge_chains = FALSE)
r_eff_seitd <- relative_eff(exp(log_lik_seitd), cores = 2)

log_lik_seeiittd <- extract_log_lik(rstan_seeiittd, merge_chains = FALSE)
r_eff_seeiittd <- relative_eff(exp(log_lik_seeiittd), cores = 2)

loo_seitd <- loo(log_lik_seitd, r_eff = r_eff_seitd, cores = 2)
loo_seeiittd <- loo(log_lik_seeiittd, r_eff = r_eff_seeiittd, cores = 2)

# Compare Models
comp <- loo_compare(loo_seitd, loo_seeiittd)
