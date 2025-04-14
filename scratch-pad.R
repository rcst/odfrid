library(data.table)
sum(softmax(1:10))

S <- 15
S*(S-1)/2
x <- c(seq(S-1, 0), seq(0, S-1))
w <- load(x)
length(w)
ll <- rod(x)
x_check <- as.numeric(routing_matrix(S) %*% ll[["y"]])
wrongs <- which(x!=x_check, arr.ind = TRUE)
wrongs
x[wrongs]
x_check[wrongs]
ll[["z"]] |> as.vector()
ll[["y"]] |> as.vector()
ll[["pi"]] |> as.vector()
ll[["q"]] |> as.vector()

x
x_check
dg <- as.data.table(ll[["dg"]])
dg[j == 2]


sum(((a <- runif(1000000)) >=0.99999)

a[which(a>=0.999999)]


# what is faster? (int)runif(...) or sample.int(1)
library(microbenchmark)
microbenchmark(as.integer(runif(1, 0, 6)), 
               sample.int(5, size = 1, replace = TRUE)-1, 
               times = 10000)

microbenchmark(rod(c(seq(20, 0), seq(0, 20))),
               rod(c(seq(40, 0), seq(0, 40))),
               rod(c(seq(75, 0), seq(0, 75))), unit = "milliseconds")









# random alighting selection
hist(rbinom(1000, 5, 0.5))
hist(sample.int(6, size = 100000, replace = TRUE)-1)
hist(as.integer(runif(n = 100000, min = 0, max = 6)))

class(my_choose(1:10, 1:10))

library(rim)
library(data.table)

maxima.start(restart = TRUE)
maxima.get("eq: j = S * (i-1) - i * (i+1) / 2 + k")
# maxima.get("eqB: idx = k - 1")
# maxima.get("assume(i>=1, i<=S)")
# maxima.get("assume(k>=i+1, k<=S)")
maxima.eval("rhs(solve([eq], k)[1])", envir = list(i = 3, j = 4))
# maxima.get("assume(equal(j,i + 1))")
# maxima.get("eqn: subst([eq], j=i+1)")
# maxima.get("solve(eq, j)")
# maxima.get("subst([S = 4, i = [1, 2, 3, 4]], y)")
maxima.eval(fk, envir = list(i = 1:5, j = 1:5))

# from stackoverflow: https://math.stackexchange.com/questions/646117/how-to-find-a-function-mapping-matrix-indices

ij_to_k <- function(i, j, n) {
  (2*n*i - i^2 + 2 * j - 3*i - 2) / 2
}

S <- 6
M <- S*(S-1)/2
dt <- data.table(i = c(rep(0, S-1),
                       rep(1, S-2),
                       rep(2, S-3),
                       rep(3, S-4),
                       rep(4, S-5)),
                 j = c(1:5, 2:5, 3:5, 4:5, 5))
dt[, k := ij_to_k(i, j, S)]
