library('cmdstanr')

printf <- function(msg, ...) cat(sprintf(msg, ...))

df <- read.csv('../data/brasil-covid-tweets.csv')
item <- df$item_id
rater <- df$annotator_id
label <- df$label_id

N <- dim(df)[1]
I <- max(item)
J <- max(rater)
K <- max(label)

# initializes pi at MLE and theta at prior
init_fun <- function(chain_id) {
    pi_tbl <- table(label) / sum(table(label))
    pi <- rep(0, K)
    for (k in 1:K) pi[k] <- pi_tbl[k]
    theta <- array(NA, c(J, K, K))
    phi <- array(NA, c(K, K))
    for (k1 in 1:K) {
        for (k2 in 1:K) {
            phi[k1, k2] <- 1 / (2 * K)
            if (k1 == k2) phi[k1, k2] <- 1 - (K - 1) / (2 * K);
        }
    }
    for (j in 1:J)
        theta[j, , ] <- phi
    list(pi = pi, theta = theta)
}

tweet_data <- list(N = N, I = I, J = J, K = K, item = item, rater = rater,
                   label = label)

printf("total labels (N) = %d;  items (I) = %d;  annotators (J) = %d;  categories (K) = %d\n",
       N, I, J, K)

model <- cmdstan_model('../stan/dawid-skene.stan')
fit <- model$sample(data = tweet_data,
                    init = init_fun,
                    chains = 4, parallel_chains = 4)

