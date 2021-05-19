library(rstan)
library(dplyr)
library(cmdstanr)
library(tibble)
library(ggplot2)

# If necessary you can load with
files <- list.files(here::here("SEIR-model", "results", "deaths_matrix_exp", "daily_seeiittd"), full.names = TRUE)

# Results
stanfit <- read_stan_csv(files)
deaths <- as_cmdstan_fit(files)

# Predicted Deaths
pred_deaths <- deaths$summary("pred_deaths")

# R_t
r_t <- deaths$summary("effective_reproduction_number")

r_t %>% 
  summarise(mean_mean = mean(mean),
            mean_median = mean(median),
            median_mean = median(mean),
            median_median = median(median))

# Omega
omega <- deaths$summary("omega")

# Real data
br <- readRDS(here::here("SEIR-model/", "data", "brazil_nation.rds"))
real_deaths <- br %>%
  pull(new_deaths) %>%
  enframe(name = "day", value = "real_deaths")

# MAE Real vs Predicted
pred_deaths %>% 
  bind_cols(real_deaths) %>% 
  mutate(
    MAE_median = abs(median - real_deaths),
    MAE_mean = abs(mean - real_deaths)) %>% 
  summarise(
    MAE_median = mean(MAE_median),
    MAE_mean = mean(MAE_mean))

# Figure Real vs Predicted
pred_deaths %>% 
  bind_cols(real_deaths) %>%
  ggplot(aes(x = day)) +
  geom_line(aes(y = median, color = "predicted"), alpha = 0.7) +
  geom_line(aes(y = real_deaths, color = "real"), alpha = 0.7) +
  geom_ribbon(aes(ymin = q5, ymax = q95), color = "grey", alpha = 0.3) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(
    title = "Daily Deaths - Predicted vs Real",
    subtitle = "Beta Regularization = 0.10",
    caption = "from 2020-02-25 to 2021-05-04",
    y = NULL,
    x = "days"
  ) +
  scale_color_brewer(palette = "Set1")

ggsave(here::here("SEIR-model", "images", "daily_deaths_betareg_010.png"), device = "png",
       dpi = 300, width = 6, height = 4)


# Plotting Stuff
stan_dens(stanfit, pars = "omega", separate_chains = TRUE)
