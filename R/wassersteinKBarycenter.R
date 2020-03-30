#' wassersteinKBarycenter
#'
#'  A wrapper for calculating K-barycenters of multivariate normal distributions with
#'   the 2-Wasserstein distance
#'
#' @param i A dummy variable ment for use with apply.
#' @param k Number k of elements in the k-barycenter.
#' @param alpha Level of trimming.
#' @param initialization Type of initialization in c('rnd', 'plus-plus'). 'rnd' makes
#'  the common random initilaization while 'plus-plus' initializes in a similar fashion to k-means++.
#' @param pooled.clusters List of multivariate normals for which the trimmed k-barycenter should be performed.
#'
#' @return A list with entries:
#' \describe{
#'   \item{wasserstein.var}{A double representing the Wasserstein variation.}
#'   \item{wasserstein.k.barycenter}{List with three elements. Variacion_wasser is Waserstein variation.
#'    Baricentro is a list of k elements, each of which is a member of the k-barycenter.
#'    Cluster is the assignation to each barycenter of the original entries.}
#' }
#'
#' @examples
#' normals <- list(list(mean = c(1, 1), cov = diag(2, 2)), list(mean = c(1, 1), cov = diag(1, 2)),
#' list(mean = c(3, 3), cov = diag(1, 2)))
#' wkb <- wassersteinKBarycenter(1, 2, 0, 'rnd', normals)
#' print(wkb$wasserstein.var)
#' print(wkb$wasserstein.k.barycente)
#'
#'@references E del Barrio, H Inouzhe, JM Loubes, C Matran and A Mayo-Iscar. (2019)
#' optimalFlow: Optimal-transport approach to flow cytometry gating and population matching.
#'  arXiv:1907.08006
#'
#' @noRd
wassersteinKBarycenter <- function(i = 1, k, alpha = 0, initialization = "rnd", pooled.clusters) {
    bar <- tryCatch(optimalFlow::trimmedKBarycenter(k = k, alpha0 = alpha, type.ini = initialization, pooled.clusters), error = function(x) list(bar = "error", 
        variacion_wasser = Inf))
    variation <- bar$variacion_wasser
    return(list(wasserstein.var = variation, wasserstein.k.barycenter = bar))
}
