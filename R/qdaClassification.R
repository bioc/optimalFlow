#' qdaClassification
#'
#' Gives quadratic discriminant scores to the points in data for a multivariate normal.
#'
#' @param normal A list with arguments mean, covaruance and weight.
#' @param data Data frame or matrix on which to perform qda.
#'
#' @return A score for each point.
#'
#' @examples
#' data.qda = cbind(rnorm(50), rnorm(50))
#' exp(qdaClassification(list(mean = c(0,0), cov = diag(1,2), weight = 1), data.qda))
#'
#' @export
#'
qdaClassification <- function(normal, data) {
    data <- as.matrix(data)
    media <- matrix(normal$mean, nrow = 1, ncol = length(normal$mean))
    inv <- tryCatch(solve(normal$cov), error = function(x) NA)
    if (is.na(inv[1])) {
        return(inv)
    } else {
        weight <- normal$weight
        qda_score <- -((1/2) * rowSums((data %*% inv) * data)) + data %*% inv %*% t(media) + matrix(-(1/2) * media %*% inv %*% 
            t(media) - (1/2) * log(1/det(inv)) + log(weight), nrow = dim(data)[1], ncol = 1)
        
    }
}
