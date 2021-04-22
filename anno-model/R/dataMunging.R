d <- read.csv("../data/long_format.csv")
d1 <- d[d$annotator_id==1,]
d2 <- d[d$annotator_id==2,]
d3 <- d[d$annotator_id==3,]


collap_d2 <- d2
collap_d2[,3] <- ifelse(collap_d2[,3]>2,1,collap_d2[,3])

#write.csv(collap_d2,file = "../data/2d.csv", row.names = FALSE)







