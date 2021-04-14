# Based on Andre's code: https://github.com/andrelmfsantos/Brazil_IO/blob/main/Brasil_IO_API.Rmd
# API doc can be found here: https://github.com/turicas/covid19-br/blob/master/api.md

library(dplyr)
library(lubridate)

download_brasilio_table <- function(dataset, table_name){
  url <- sprintf("https://data.brasil.io/dataset/%s/%s.csv.gz", dataset, table_name)
  tmp <- tempfile()
  download.file(url, tmp)
  response <- read.csv(gzfile(tmp), encoding = "UTF-8")
  unlink(tmp)
  return(response)
}

# table name is one of "caso", "caso_full", "obito_cartorio"
# we are choosing "caso_full" for the full data

data <- download_brasilio_table("covid19", "caso_full") # 42mb

# Getting only 2021
y21 <- data %>%
  as_tibble %>% 
  mutate(date = as_date(date)) %>% 
  filter(date <= max(date) & date >= "2021-01-01")

# Getting only state stuff
state <- y21 %>% 
  filter(city == "") %>% 
  select(-city)

# Getting Brazil's national-level data
br <- state %>% 
  group_by(date) %>% 
  summarise(
    estimated_population_2019 = sum(estimated_population_2019),
    last_available_confirmed_per_100k_inhabitants = sum(last_available_confirmed_per_100k_inhabitants),
    last_available_deaths = sum(last_available_deaths),
    new_confirmed = sum(new_confirmed),
    new_deaths = sum(new_deaths)
  )

# save file
saveRDS(state, here::here("SEIR-model", "data", "state_2021.rds"))
saveRDS(br, here::here("SEIR-model", "data", "brazil_2021.rds"))
