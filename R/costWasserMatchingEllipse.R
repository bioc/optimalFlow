#' costWasserMatchinEllipse
#'
#' Calculates a similarity distance based on the 2-Wassertein distance between mixtures of multivariate normal distributions.
#'
#' @param test.cytometry A clusetering represented as a list of clusters. Each cluster is a list with elements mean, cov, weight and type.
#' @param training.cytometries A list of clusterings with the same format as test.cytometry.
#' @param equal.weights If True, weights assigned to every cluster in a partion are uniform (1/number of clusters) when calculating the similarity distance. If False, weights assigned to clusters are the proportions of points in every cluster compared to the total amount of points in the partition.
#'
#' @return A vector representing the similarity distance between test.cytometry and the elements in training.cytometries.
#'
#' @examples
#' partition1 <- list(list(mean = c(1, 1), cov = diag(1, 2), weight = 0.5, type = '1'),
#'                   list(mean = c(-1, -1), cov = diag(1, 2), weight = 0.5, type = '2'))
#' partition2 <- list(list(list(mean = c(1, -1), cov = diag(1, 2),
#'                             weight = 0.5, type = '1'), list(mean = c(-1, 1), cov = diag(1, 2), weight = 0.5, type = '2')))
#' costWasserMatchingEllipse(partition1, partition2)
#'
#' @references E del Barrio, H Inouzhe, JM Loubes, C Matran and A Mayo-Iscar. (2019) optimalFlow: Optimal-transport approach to flow cytometry gating and population matching. arXiv:1907.08006
#'
#' @export
#'
costWasserMatchingEllipse <- function(test.cytometry, training.cytometries, equal.weights = FALSE) {
    dim.test.cytometries <- length(training.cytometries)
    cost.vector <- array(0, dim = c(1, dim.test.cytometries))
    training.cytometries_elipses <- list()
    for (i in 1:(dim.test.cytometries)) {
        training.cytometry <- training.cytometries[[i]]
        names.a <- unlist(lapply(test.cytometry, function(x) x$type))
        names.b <- unlist(lapply(training.cytometry, function(x) x$type))
        weights.a <- unlist(lapply(test.cytometry, function(x) x$weight))
        if (abs(1 - sum(weights.a)) > 10^(-5)) {
            message("Weights do not add to 1. A normalization will be applied.")
            weights.a <- weights.a/sum(weights.a)
        }
        weights.b <- unlist(lapply(training.cytometry, function(x) x$weight))
        if (abs(1 - sum(weights.b)) > 10^(-5)) {
            message("Weights do not add to 1. A normalization will be applied.")
            weights.b <- weights.b/sum(weights.b)
        }

        if (equal.weights) {
            a <- rep(1/length(names.a), length(names.a))
            b <- rep(1/length(names.b), length(names.b))
            names(a) <- names.a
            names(b) <- names.b
        } else {
            a <- weights.a
            b <- weights.b
            names(a) <- names.a
            names(b) <- names.b
        }

        nn <- length(a)
        mm <- length(b)
        cost.ab <- array(dim = c(nn, mm))

        naive.transport.cost <- 0
        for (k in 1:nn) {
            for (l in 1:mm) {
                cost.ab[k, l] <- optimalFlow::w2dist(test.cytometry[[k]], training.cytometry[[l]])
                naive.transport.cost <- naive.transport.cost + cost.ab[k, l] * a[k] * b[l]
            }
        }
        tranport.plan <- transport::transport(a, b, cost.ab)
        transport.cost <- 0
        for (k in 1:dim(tranport.plan)[1]) {
            transport.cost <- transport.cost + cost.ab[tranport.plan[k, 1], tranport.plan[k, 2]] * tranport.plan[k, 3]
        }
        cost.vector[1, i] <- transport.cost/naive.transport.cost
    }
    return(cost.vector)
}
