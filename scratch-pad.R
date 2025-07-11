i_to_id(i = 0, S = 6)
i_to_id(i = 0, N = 3, S = 6)

library(data.table)
sum(softmax(1:10))


S <- 40
N <- 10
D <- S*(S-1)/2
x <- rep(x = c(seq(S-1, 0), seq(0, S-1)), times = N) |> matrix(ncol = N)
# checking for rare negative OD
load(x)
ll <- rod(x)
ll[["y"]]
ll[["z"]]
ll[["q"]]
x_check <- as.numeric(routing_matrix(S) %*% ll[["y"]]) |> matrix(ncol = N)



all(x == x_check)

for(i in 1:10) {
  ll <- rod_slow(x)
  stopifnot(all(ll[["y"]] >=0))
  x_check <- as.numeric(routing_matrix(S) %*% ll[["y"]])
  stopifnot(all(x == x_check))
}


# what is faster? (int)runif(...) or sample.int(1)
library(microbenchmark)
microbenchmark(as.integer(runif(1, 0, 6)), 
               sample.int(5, size = 1, replace = TRUE)-1, 
               times = 10000)

microbenchmark(rod(c(seq(20, 0), seq(0, 20))),
               rod(c(seq(40, 0), seq(0, 40))),
               rod(c(seq(75, 0), seq(0, 75))),
               rod(c(seq(100,0), seq(0, 100))), unit = "seconds")









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


# random intergers that are lower than max
v <- 20
y <- ivrunif(min = 0, max = c(5, 10, 3, 7, 14), sum =0) 
u <- -1.0 * log(y)
p <- u/sum(u) * v
sum(p)

class(uniformSimplexSample(N = 20, C = 5))

z <- c(5, 3, 7, 8, 8, 8, 10)
x <- c(4, 3, 6, 9, 3, 7, 12)
a <- adjustWithCaps(x = x, z = z)
sum(x) == sum(a)

S <- 5
S*(S-1)/2.0
M <- (S*(S-1)/2.0) - 1
N <- 3
D <- (1.0 + sqrt(1.0 + 8.0 * (M+1))) / 2.0
S == D

G <- matrix(data = rnorm(n = M*N), nrow = M, ncol = N)
matrix_softmax(G, 1.0)

microbenchmark(arma = log_choose_vector(n = 1:1000, k = rep(3, 1000)),
            R = lchoose(n = 1:1000, k = rep(3, 1000)))

log_choose(N = matrix(data = rep(1:1000, 3), ncol = 3), K = matrix(data = rep(rep(3, 1000), 3), ncol = 3)) |> exp()

odformsub2ind(6, 4, 2)
test(6, 2, 3)
