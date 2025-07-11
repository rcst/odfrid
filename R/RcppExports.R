# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

adjustWithCaps <- function(x, z) {
    .Call(`_odfrid_adjustWithCaps`, x, z)
}

#' @param k covariance matrix for column psi_d
#' @param psi mapping-factor matrix
#' @param d the index of column of matrix psi to be sampled
NULL

#' load - No. passengers on the bus immediatly after each stops
#'
#' @param x column vector of boardings and alightings
#' @return vector of passengers loadings immediatly after each stops
load <- function(x) {
    .Call(`_odfrid_load`, x)
}

#' routing_matrix - construct routing matrix from no. stops (S)
#'
#' @param S length-one integer denoting the number of stops
#' @return An (armadillo) integer matrix
routing_matrix <- function(s) {
    .Call(`_odfrid_routing_matrix`, s)
}

#' ztoy - capped uniform simplex sampling 
#'
#' Helper function that generates an integer vector with constant total sum and varying caps
#'
#' @param z numeric vector of no. passengers approaching a single specific stop
#' @param v double no. alighters at that specific stop
#' @return An integer vector of same length as z whos values are all smaller-or-equal to z and it's sum is equal to v
ztoy <- function(z, v) {
    .Call(`_odfrid_ztoy`, z, v)
}

#' rod - Conditional Sampling of OD vectors
#' 
#' @param x a integer matrix whoes columns each contain boarding and alighting counts of 1 bus
#' journey
#' @return A named list of containing (1) the sampled OD vector (named y), (2) a corresponging vector (named z)
#' the log probability density from Markov chain transition probabilities (named lq)
rod <- function(x) {
    .Call(`_odfrid_rod`, x)
}

roundWithPreservedSum <- function(fn) {
    .Call(`_odfrid_roundWithPreservedSum`, fn)
}

uniformSimplexSample <- function(N, C = 1.0, A = 0.0) {
    .Call(`_odfrid_uniformSimplexSample`, N, C, A)
}

