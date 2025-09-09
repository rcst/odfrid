#' generate_fake_pc - Generate Fake Passenger Counts
#'
#' @param N the number of trips
#' @param S the number of stops
#' @param boarders The average number of boarder per stop.
#' @return A matrix where each column is an OD vector
#' @export
generate_fake_pc <- function(N = 100, t_step = 300, S = 6, dest, boarders = 20) {
  M <- S * (S - 1) / 2

  if(missing(dest))
     dest <- sample(x = 1:S, size = 1) 

  # index-tracking matrix
  K <- matrix(nrow = S, ncol = S)
  K[lower.tri(K)] <- 1:M
  K <- t(K)

  i <- function(j) K[,j][1:(j-1)]
  j <- function(i) K[i,][(i+1):S]

  lbd <- matrix(data = 0.1, nrow = M, ncol = N)

  # set all alighting on last stops to 1
  lbd[M,] <- 1.0
  pA <- 0.9
  pB <- 0.1
  for(t in 1:N) {
    pA_t <- (pB - pA) * (t-1)/N + pA
    pB_t <- (pA - pB) * (t-1)/N + pB
    # lbd[i(dest[1]),t] <- pA_t
    lbd[i(dest),t] <- pB_t

    # normalize to 1
    for(k in 1:(S-2)) 
	    lbd[j(k), t] <- lbd[j(k),t] / sum(lbd[j(k),t])
  }

  # generate OD vector
  y <- matrix(data = NA, nrow = M, ncol = N)
  for(t in 1:N) {
    for(k in 1:(S-1)) {
      # brds <- as.integer(rnorm(n = 1, mean = sqrt(boarders), sd = sqrt(5))^2)
      y[j(k),t] <- rmultinom(n = 1, size = boarders, prob = lbd[j(k),t])
    }
  }

  # calculate passenger counts
  x <- routing_matrix(S) %*% y

  return(list(S = S, 
              destinations = dest, 
              t =  (0:(N-1)) * t_step,
              x = x, 
              y = y, 
              lambda = lbd))
}

#' @export
j2i <- function(S) {
  M <- S * (S-1) / 2
  K <- matrix(nrow = S, ncol = S)
  K[lower.tri(K)] <- 1:M
  K <- t(K)

  K[i,][(i+1):S]
}

#' @export
i2j <- function(S) {
  M <- S * (S-1) / 2
  K <- matrix(nrow = S, ncol = S)
  K[lower.tri(K)] <- 1:M
  K <- t(K)

  K[,j][1:(j-1)]
}

#' @export
idx2ij <- Vectorize(FUN = function(S, idx) {
			    M <- S * (S-1) / 2
			    K <- matrix(nrow = S, ncol = S)
			    K[lower.tri(K)] <- 1:M
			    K <- t(K)
			    which(K == idx, arr.ind = TRUE) 
}, vectorize.args = "idx")
