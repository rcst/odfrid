library(data.table)
sum(softmax(1:10))

S <- 75 
S*(S-1)/2
x <- c(seq(S-1, 0), seq(0, S-1))
# checking for rare negative OD
ll <- rod(x)
for(i in 1:10) {
  ll <- rod(x)
  stopifnot(all(ll[["y"]] >=0))
  x_check <- as.numeric(routing_matrix(S) %*% ll[["y"]])
  stopifnot(all(x == x_check))
}

library(microbenchmark)
microbenchmark(rod(c(seq(20, 0), seq(0, 20))),
               rod(c(seq(40, 0), seq(0, 40))),
               rod(c(seq(75, 0), seq(0, 75))),
               rod(c(seq(100,0), seq(0, 100))), unit = "seconds")
