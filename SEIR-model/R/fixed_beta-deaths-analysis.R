library(dplyr)
library(cmdstanr)
library(tibble)
library(ggplot2)

# If necessary you can load with
deaths_seitd <- list.files(here::here("SEIR-model", "results", "deaths_fixed_beta", "seitd"), full.names = TRUE) %>% 
  as_cmdstan_fit()
deaths_seeiittd <- list.files(here::here("SEIR-model", "results", "deaths_fixed_beta", "seeiittd"), full.names = TRUE) %>% 
  as_cmdstan_fit()

# Predicted Deaths
pred_deaths_seitd <- deaths_seitd$summary("pred_deaths")
pred_deaths_seeiittd <- deaths_seeiittd$summary("pred_deaths")

# Real data
br <- readRDS(here::here("SEIR-model/", "data", "brazil_nation.rds"))
real_deaths <- br %>%
  group_by(week = cut(date, "week")) %>%
  summarise(deaths = sum(new_deaths)) %>%
  pull(deaths) %>%
  ceiling %>%
  as.integer %>% 
  head(-1) %>% 
  enframe(name = "week", value = "real_deaths")

# MAE Real vs Predicted
pred_deaths_seitd %>% 
  bind_cols(real_deaths) %>% 
  mutate(
    MAE_median = abs(median - real_deaths),
    MAE_mean = abs(mean - real_deaths)) %>% 
  summarise(
    MAE_median = mean(MAE_median),
    MAE_mean = mean(MAE_mean))

pred_deaths_seeiittd %>% 
  bind_cols(real_deaths) %>% 
  mutate(
    MAE_median = abs(median - real_deaths),
    MAE_mean = abs(mean - real_deaths)) %>% 
  summarise(
    MAE_median = mean(MAE_median),
    MAE_mean = mean(MAE_mean))

# Figure Real vs Predicted
real_deaths %>% 
  bind_cols(pred_deaths_seitd %>% select(median_seitd = median)) %>% 
  bind_cols(pred_deaths_seeiittd %>% select(median_seeiittd = median)) %>% 
  ggplot(aes(x = week)) +
  geom_line(aes(y = median_seitd, color = "SEITD predicted"), alpha = 0.7) +
  geom_line(aes(y = median_seeiittd, color = "SEEIITTD predicted"), alpha = 0.7) +
  geom_line(aes(y = real_deaths, color = "real"), alpha = 0.7) +
  #geom_ribbon(aes(ymin = q5, ymax = q95), color = "grey", alpha = 0.3) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(
    title = "Weekly Deaths",
    subtitle = "Predicted vs Real",
    caption = "from 2020-02-25 to 2021-05-04",
    y = NULL,
    x = "weeks"
  ) +
  scale_color_brewer(palette = "Set1")

ggsave(here::here("SEIR-model", "images", "fixed_beta_weekly_deaths.png"), device = "png",
       dpi = 300, width = 6, height = 4)
