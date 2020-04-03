#' labelTransferEllipse
#'
#' Label transfer between a test partition and a training partitions viewed as a mixture of gaussians.
#'
#' @param i A dummy variable, should be any integral. Ment for use with lapply.
#' @param test.cytometry.ellipses A test clustering viewed as a mixture of multivariate normal distributions.
#' @param training.cytometries.barycenter A training partition viewed as a mixture of multivariate normal distributions.
#' @param equal.weights If True, weights assigned to every cluster in a partion are uniform (1/number of clusters) when calculating the similarity distance. If False, weights assigned to clusters are the proportions of points in every cluster compared to the total amount of points in the partition.
#'
#' @return A fuzzy relabeling consistent of a transportation plan.
#'
#' @examples
#' partition1 <- list(list(mean = c(1, 1), cov = diag(1, 2), weight = 0.5, type = '1'),
#'                   list(mean = c(-1, -1), cov = diag(1, 2), weight = 0.5, type = '2'))
#' partition2 <- list(list(mean = c(1, 1), cov = diag(1, 2), weight = 0.5, type = 'a'),
#'                   list(mean = c(-1, -1), cov = diag(1, 2), weight = 0.5, type = 'b'))
#' labelTransferEllipse(1, partition2, partition1)
#' @references E del Barrio, H Inouzhe, JM Loubes, C Matran and A Mayo-Iscar. (2019) optimalFlow: Optimal-transport approach to flow cytometry gating and population matching. arXiv:1907.08006
#'
#' @export
#'
labelTransferEllipse <- function(i, test.cytometry.ellipses, training.cytometries.barycenter, equal.weights = FALSE) {
    N <- length(training.cytometries.barycenter)
    M <- length(test.cytometry.ellipses)
    cost.matrix <- array(dim = c(N, M))
    for (j in 1:N) {
        for (i in 1:M) {
            cost.matrix[j, i] <- w2dist(list(mean = test.cytometry.ellipses[[i]]$mean, cov = test.cytometry.ellipses[[i]]$cov),
                list(mean = training.cytometries.barycenter[[j]]$mean, cov = training.cytometries.barycenter[[j]]$cov))
        }
    }
    names.a <- unlist(lapply(training.cytometries.barycenter, function(x) x$type))
    names.b <- unlist(lapply(test.cytometry.ellipses, function(x) x$type))

    if (equal.weights) {
        A <- matrix(1/N, nrow = 1, ncol = N)
        B <- matrix(1/M, nrow = 1, ncol = M)
        names(A) <- names.a
        names(B) <- names.b
    } else {
        A <- unlist(lapply(training.cytometries.barycenter, function(x) x$weight))
        A <- A/sum(A)
        names(A) <- names.a
        B <- unlist(lapply(test.cytometry.ellipses, function(x) x$weight))
        names(B) <- names.b
    }
    optimal.transport.form.A.to.B <- transport::transport(a = A, b = B, costm = cost.matrix)
    optimal.transport.form.A.to.B$to <- names(B)[optimal.transport.form.A.to.B$to]
    optimal.transport.form.A.to.B$from <- names(A)[optimal.transport.form.A.to.B$from]
    return(optimal.transport.form.A.to.B)
}
