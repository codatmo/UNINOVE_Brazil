library(cmdstanr)
library(rstan)

d <- read.csv("../data/long_format.csv")

d1 <- d[d$annotator_id==1,]
d2 <- d[d$annotator_id==2,]
d3 <- d[d$annotator_id==3,]


cd2 <- d2
cd2[,3] <- ifelse(cd2[,3]>2,1,cd2[,3]) # collapse 3-8 to 2

# split on white space
cd2$words <- strsplit(as.character(cd2$text)," ")


# make data factors
# fit with logistic regression


