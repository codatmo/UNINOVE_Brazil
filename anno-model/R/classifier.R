library(cmdstanr)
# Description of this code is at: 
# https://github.com/codatmo/UNINOVE_Sao_Paulo/tree/main/anno-model/TweetClassifier.html

d <- read.csv("../data/long_format.csv") # Not included with code due to IP restrictions

d1 <- d[d$annotator_id==1,]
d2 <- d[d$annotator_id==2,]
d3 <- d[d$annotator_id==3,]

cd2 <- d2
cd2$label_id <- ifelse(cd2$label_id > 2, 1, cd2$label_id) # collapse annotations 3-8 to 2
cd2$label_id <- cd2$label_id - 1 # map 2 to 1, 1 to 0
# split on white space
cd2$words <- strsplit(as.character(cd2$text),"\\s+", perl = TRUE)

words <- unlist(cd2$words)
unique_words <- unique(words)
vocab <- new.env(size = length(unique_words))
counts <- new.env(size = length(unique_words))
index <- 0L

# ugly but base R, there is a hash package 
for (i in 1:length(words)) {
  if (!exists(words[i], envir = vocab, inherits = FALSE)) {
    index <- index + 1
    assign(words[i], index, envir = vocab)
    assign(words[i], 1, envir = counts)
  }
  else {
    increment <- get(words[i], envir = counts, inherits = FALSE)
    assign(words[i], increment + 1, envir = counts)
  }
}

# return vector of length vocab with words in tweet set to 1, others 0
words2binary_vector <- function(words, vocab_lookup) {
  words_unlisted <- unlist(words)
  sent_vect <- rep.int(0,length(vocab_lookup))
  for (i in 1:length(words_unlisted)) {
     if (exists(words_unlisted[i], envir = vocab_lookup, inherits = FALSE)) {
       sent_vect[get(words_unlisted[i], envir = vocab_lookup, inherits = FALSE)] <- 1
     }
  }
  return(sent_vect)
}

binary_pred_matrix <- matrix(nrow = nrow(cd2), 
       ncol = index)

for (i in 1:nrow(cd2)) {
  binary_pred_matrix[i,] <- words2binary_vector(cd2[i,]$words, vocab_lookup = vocab)
}

N_data_count <- nrow(binary_pred_matrix)
held_out <- 100
stan_data <- list(N_data_count = N_data_count - held_out + 1, 
                  M_feature_count = ncol(binary_pred_matrix),
                  H_held_out = held_out,
                  x_features = binary_pred_matrix[held_out:N_data_count,], 
                  y_categories = cd2[held_out:N_data_count,]$label_id,
                  x_held_out = binary_pred_matrix[1:held_out,],
                  y_held_out = cd2[1:held_out,]$label_id)

model <- cmdstan_model("../stan/classifier.stan")

fit <- model$sample(data = stan_data, seed = 22, chains = 4, 
                    iter_warmup = 1000, iter_sampling = 1000)

# evaluating on training data, not a good idea
print(fit$draws(c('precision', 'recall')))

fit$save_object("held_out1.RDS")

library(shinystan)
library(rstan)

stanfit <- rstan::read_stan_csv(fit$output_files())
launch_shinystan(stanfit)


# fit_mle <- model$optimize(data = stan_data, seed = 22)
