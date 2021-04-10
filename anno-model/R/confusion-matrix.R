df <- read.csv('../data/brasil-covid-tweets.csv')
item <- df$item_id
rater <- df$annotator_id
label <- df$label_id

I <- max(item)
K <- max(label)

# calculate confusion matrix  
conf <- matrix(0, K, K)
for (i in 1:I) {
    labs <- label[item == i]
    L <- length(labs)
    for (j1 in 1:L) {
        k1 <- labs[j1]
        for (j2 in 1:L) {
            k2 <- labs[j2]
            conf[k1, k2] <- conf[k1, k2] + 1
        }
    }
}
cat("\nCo-annotation matrix\n")
print(conf)

# agreement by category --- have to print this on the outside
agree <- array(0, K)
for (k in 1:K)
    agree[k] <- conf[k, k] / sum(conf[k , ])

cat("\nAgreement by category\n")
print(agree, digits = 2)
