library(cmdstanr)
#library(rstan)
#library(rstanarm)

d <- read.csv("../data/long_format.csv")

d1 <- d[d$annotator_id==1,]
d2 <- d[d$annotator_id==2,]
d3 <- d[d$annotator_id==3,]


cd2 <- d2
cd2[,3] <- ifelse(cd2[,3]>2,1,cd2[,3]) # collapse 3-8 to 2

# split on white space
cd2$words <- strsplit(as.character(cd2$text),"\\s+", perl = TRUE)

vocab <- new.env(hash = TRUE)
counts <- new.env(hash = TRUE)
index <- 0L
words <- unlist(cd2$words)

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

words2coefs <- function(words, index_lookup) {
  words_unlisted <- unlist(words)
  coefs <- rep.int(0,length(index_lookup))
  for (i in 1:length(words_unlisted)) {
    coefs[get(words_unlisted[i], envir = index_lookup, inherits = FALSE)] <- 1
  }
  return(coefs)
}

coefs <- matrix(nrow = nrow(cd2), 
       ncol = length(vocab))

for (i in 1:nrow(cd2)) {
  coefs[i,] <- words2coefs(cd2[i,]$words, index_lookup = vocab)
}

stan_data <- list(N_data_count = nrow(coefs), 
                  M_feature_count = ncol(coefs),
                  x_features = coefs, 
                  y_categories = cd2$label_id)

model <- cmdstan_model("classifier.stan")

fit <- model$sample(data = stan_data, seed = 22, chains = 1)


