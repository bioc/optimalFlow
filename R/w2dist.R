#' w2dist
#'
#' The 2-Wasserstein distance between two multivariate normal distributions
#'
#' @param P A multivariate normal distribution given as a list with arguments mean and cov.
#' @param Q A multivariate normal distribution given as a list with arguments mean and cov.
#'
#' @return A double giving the 2-Wasserstein distance between the two distributions.
#'
#' @examples
#' P <- list(mean = c(1, 1), cov = diag(1, 2))
#' Q <- list(mean = c(0, 0), cov = 1.1*diag(1, 2))
#' Q <- list(mean = c(0, 0), cov = 1.1*diag(1, 2))
#' w2dist(P, Q)
#'
#' @export
#'
w2dist <- function(P, Q) {
    sqrt(abs(optimalFlow:::distGaussianMean(P, Q) + optimalFlow:::distGaussianCov(P, Q)))
}
