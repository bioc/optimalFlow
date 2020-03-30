#' estimCovCellGeneral
#'
#' Estimation of mean and covariance for a label in a partition.
#'
#' @param cell Labell of the clsuter of interest.
#' @param cytometry Data of the partition, without labels.
#' @param labels Labels of the partition.
#' @param type How to estimate covariance matrices of a cluster. 'standard' is for using cov(), while 'robust' is for using robustbase::covMcd.
#' @param alpha Only when type = 'robust'. Indicates the value of alpha in robustbase::covMcd.
#'
#' @return A list containing:
#' \describe{
#'  \item{mean}{Mean of the cluster.}
#'  \item{cov}{Covariance of the cluster.}
#'  \item{weight}{Weight associated to the cluster.}
#'  \item{type}{Type of the cluster.}
#' }
#'
#' @examples
#' estimCovCellGeneral('Basophils', Cytometry1[,1:10], Cytometry1[,11])
#'
#' @importFrom stats cov
#'
#' @export
#'
estimCovCellGeneral <- function(cell, cytometry, labels, type = "standard", alpha = 0.85) {
    indices <- which(labels == cell)
    ddd <- dim(cytometry)
    size <- length(indices)
    weight <- size/ddd[1]
    cytometry <- cytometry[indices, 1:(ddd[2])]
    if (size >= 2) {
        if (type == "standard") {
            return(list(mean = colMeans(cytometry), cov = cov(cytometry), weight = weight, type = cell))
        } else {
            if (type == "robust" & alpha * size > (ddd[2] + 1)) {
                rob_est <- tryCatch(robustbase::covMcd(cytometry, alpha = alpha), error = function(x) NA)
                if (is.na(rob_est)[1]) {
                  return(NA)
                } else {
                  return(list(mean = rob_est$center, cov = rob_est$cov, weight = weight, type = cell))
                }
            } else {
                return(NA)
            }
        }
    } else {
        return(NA)
    }
}