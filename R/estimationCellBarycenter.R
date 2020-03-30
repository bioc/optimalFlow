#' estimationCellBarycenter
#'
#' Estimates a Wasserstein barycenter for a cluster type using a collection of partitions.
#'
#' @param cell Name of the cluster of interest.
#' @param cytometries List of clusterings.
#'
#' @return A list representing the (1-)barycenter:
#' \describe{
#'  \item{mean}{Mean of the barycenter.}
#'  \item{cov}{Covariance of the barycenter.}
#'  \item{weight}{Weight associated to the barycenter.}
#'  \item{type}{Type of the cluster.}
#' }
#'
#' @examples
#' partition1 <- list(list(mean = c(1, 1), cov = diag(1, 2), weight = 0.5, type = '1'),
#'                   list(mean = c(-1, -1), cov = diag(1, 2), weight = 0.5, type = '2'))
#' partition2 <- list(list(mean = c(1, -1), cov = diag(1, 2), weight = 0.5, type = '1'),
#'                   list(mean = c(-1, 1), cov = diag(1, 2), weight = 0.5, type = '2'))
#' cytometries <- list(partition1, partition2)
#' estimationCellBarycenter('1',cytometries)
#'
#' @export
#'
estimationCellBarycenter = function(cell, cytometries) {
    inexes.cell <- lapply(cytometries, function(x) which(unlist(lapply(x, function(y) y$type)) == cell))
    cell.representatives <- list()
    j <- 0
    for (i in 1:length(inexes.cell)) {
        if (length(inexes.cell[[i]]) == 0) {
            next
        } else {
            j <- j + 1
            cell.representatives[[j]] <- cytometries[[i]][[inexes.cell[[i]]]]
        }
    }
    if (length(cell.representatives) > 0) {
        weight.cell <- unlist(lapply(cell.representatives, function(x) x$weight))
        cov.cell <- lapply(cell.representatives, function(x) x$cov)
        cov.barycenter.cell <- optimalFlow:::GaussianBarycenters(cov.cell, weight.cell/sum(weight.cell))$Barycenter
        mean.cell <- lapply(cell.representatives, function(x) x$mean)
        mean.barycenter.cell <- colSums(matrix(unlist(lapply(1:length(cell.representatives), function(i) mean.cell[[i]] * weight.cell[i]/sum(weight.cell))),
            nrow = length(cell.representatives), byrow = TRUE))
        weight.barycenter.cell = mean(weight.cell)
        return(list(mean = mean.barycenter.cell, cov = cov.barycenter.cell, weight = weight.barycenter.cell, type = cell))
    }
}
