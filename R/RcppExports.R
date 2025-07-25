# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

broadcast_test_1 <- function() {
    .Call(`_odfrid_broadcast_test_1`)
}

broadcast_test_2 <- function() {
    .Call(`_odfrid_broadcast_test_2`)
}

broadcast_test <- function() {
    .Call(`_odfrid_broadcast_test`)
}

#' routing_matrix - construct routing matrix from no. stops (S)
#'
#' @param S length-one integer denoting the number of stops
#' @return An (armadillo) integer matrix
routing_matrix <- function(s) {
    .Call(`_odfrid_routing_matrix`, s)
}

model_sample <- function(ax, dep_time, sample, warmup, D, print_n = 100L) {
    .Call(`_odfrid_model_sample`, ax, dep_time, sample, warmup, D, print_n)
}

#' ztoy - capped uniform simplex sampling 
#'
#' Helper function that generates an integer vector with constant total sum and varying caps
#'
#' @param z numeric vector of no. passengers approaching a single specific stop
#' @param v double no. alighters at that specific stop
#' @return An integer vector of same length as z whos values are all smaller-or-equal to z and it's sum is equal to v
NULL

#' rod - Conditional Sampling of OD vectors
NULL

#' load - No. passengers on the bus immediatly after each stops
#'
#' @param x column vector of boardings and alightings
#' @return vector of passengers loadings immediatly after each stops
load <- function(x) {
    .Call(`_odfrid_load`, x)
}

rod <- function(x) {
    .Call(`_odfrid_rod`, x)
}

#' ess_psi - Elliptical Slice Sampling of Psi
#'
#' One call to this function wll update the global matrix Psi based on the 
#' current (global) set of parameters.
#'
#' @param k covariance matrix for column psi_d
NULL

test <- function() {
    invisible(.Call(`_odfrid_test`))
}

