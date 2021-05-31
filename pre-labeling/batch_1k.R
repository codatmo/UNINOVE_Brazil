library(googledrive)
library(dplyr, warn.conflicts = FALSE)
library(readr, warn.conflicts = FALSE)
library(tidyr, warn.conflicts = FALSE)
library(vroom)
library(fs)

# Auth if necessary
# drive_auth()

file <- drive_download(as_id("14FuBSFHX9-NyhEQhgvU-O5VdIpe23r2p"))

df <- vroom(file$local_path,
            col_names = c("index", "id", "created_at", "text", "user",
                          "place", "user_place", "country", "coordinates",
                          "undefined_col", "undefined_col2", "undefined_col3"),
            col_types = "ffTcfccc?cii",
            na = c("", "NA", "None"),
            locale = locale(date_format = "y-m-d H:M:S"))

#### Sampling from the 2.6mi BR tweets ####

# First round of 1k cross-anotation
set.seed(123)
sample_df <- df %>%
  drop_na(text) %>%
  sample_n(1e3)

# Second round of 1k cross-anotation
set.seed(124)
sample_df2 <- df %>%
  drop_na(text) %>%
  sample_n(1e3) %>%
  anti_join(sample_df)

# 100 samples to teach anotations
set.seed(125)
sample_df3 <- df %>%
  drop_na(text) %>%
  sample_n(1e2) %>%
  anti_join(sample_df) %>%
  anti_join(sample_df2)

# Third round of 1k anotation
set.seed(126)
sample_df4 <- df %>%
  drop_na(text) %>%
  sample_n(4002) %>%
  anti_join(sample_df) %>%
  anti_join(sample_df2) %>%
  anti_join(sample_df3)

# Fourth round of 1k anotation
set.seed(127)
sample_df5 <- df %>%
  drop_na(text) %>%
  sample_n(1e3) %>%
  anti_join(sample_df) %>%
  anti_join(sample_df2) %>%
  anti_join(sample_df3) %>%
  anti_join(sample_df4)

# Fifth round of 1k anotations (I need 3k for 3 labellers 1k each)
set.seed(128)
sample_df6 <- df %>%
  drop_na(text) %>%
  sample_n(3009) %>%
  anti_join(sample_df) %>%
  anti_join(sample_df2) %>%
  anti_join(sample_df3) %>%
  anti_join(sample_df4) %>%
  anti_join(sample_df5)

file_delete(file$local_path)

write_csv(sample_df6, "sample_1k.csv")
