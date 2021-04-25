library(cmdstanr)
#library(rstan)

d <- read.csv("../data/long_format.csv")

d1 <- d[d$annotator_id==1,]
d2 <- d[d$annotator_id==2,]
d3 <- d[d$annotator_id==3,]

cd2 <- d2
cd2$label_id <- ifelse(cd2$label_id > 2, 1, cd2$label_id) # collapse 3-8 to 2
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

binary_pred_vector <- matrix(nrow = nrow(cd2), 
       ncol = index)

for (i in 1:nrow(cd2)) {
  binary_pred_vector[i,] <- words2binary_vector(cd2[i,]$words, vocab_lookup = vocab)
}

stan_data <- list(N_data_count = nrow(binary_pred_vector), 
                  M_feature_count = ncol(binary_pred_vector),
                  x_features = binary_pred_vector, 
                  y_categories = cd2$label_id)

model <- cmdstan_model("../stan/classifier.stan")

#fit <- model$sample(data = stan_data, seed = 22, chains = 1)

#fit$save_object("one_chain.RDS")

library(bayesplot)

fit_mle <- model$optimize(data = stan_data, seed = 22)
