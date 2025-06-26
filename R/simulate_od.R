#' generate_od_matrix - Function to generate an O-D matrix using gravity model
#'
#' @param production
#' @param attraction
#' @param distance_matrix
#' @param total_trips
#' @param distance_decay Numeric vector of length 1 that controls the likelihood of a trip being made with increasing distance
#' @return A matrix object
#' @export
generate_od_matrix <- function(num_zones, production, attraction, distance_matrix, total_trips, distance_decay) {
  od_matrix <- matrix(0, nrow = num_zones, ncol = num_zones)
  
  for (i in 1:num_zones) {
    for (j in 1:num_zones) {
      od_matrix[i, j] <- (production[i] * attraction[j]) / (distance_matrix[i, j]^distance_decay)
    }
  }
  
  # Normalize the O-D matrix to match the total number of trips
  od_matrix <- od_matrix / sum(od_matrix) * total_trips
  return(od_matrix)
}

#' random_od_matrix - Randomly generate O-D matrices using gravity model
#'
#' @param N
#' @param stops
#' @param total_trips
#' @param distance_decay
#' @return List of O-D matrices
#' @export
random_od_matrix <- function(N = 1, stops = 6, total_trips = 200, distance_decay = 1.5) {
  # TODO
  # - generate integer values matrices
  # - add top-level list that contains
  #   - OD matrices as quadratic matrices/ upper triangular matrices
  #   - OD matrices where rows are individual OD matrices as stacked vectors and columns are individual trips
  #   - boarding and alighting vectors


  production <- runif(stops, min = 100, max = 500)  # Initial random production values
  attraction <- runif(stops, min = 100, max = 500)  # Initial random attraction values

  # can be improved
  distance_matrix <- matrix(runif(stops^2, min = 1, max = 10), nrow = stops)
  ods <- list()  # List to store O-D matrices

  for(k in 1:N) {
    ods[[k]] <- generate_od_matrix(stops, production, attraction, distance_matrix, total_trips, distance_decay)
    # Smoothly vary production and attraction values for the next time step
    production <- production + rnorm(stops, mean = 0, sd = 10)  # Small random perturbation
    attraction <- attraction + rnorm(stops, mean = 0, sd = 10)  # Small random perturbation
  }

  return(ods)
}
