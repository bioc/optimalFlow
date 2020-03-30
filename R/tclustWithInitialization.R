#' tclustWithInitialization
#'
#' A wrapper for the function tclust_H.
#'
#' @param initialization Initial solution for parameters provided by the user. Can be a matrix of data containing observations anc cluster assignations or can be a list spesifying a multivariate mixture of gaussians.
#' @param cytometry A matrix or data.frame of dimension n x p, containing the observations (row-wise).
#' @param i.sol.type Type of initial solutions in c('points', 'barycenters'). 'points' refers to a classified data matrix, while 'barycenters' to a multivariate mixture.
#' @param trimming The proportion of observations to be trimmed.
#' @param restr.fact The constant restr.fact >= 1 constrains the allowed differences among group scatters. Larger values imply larger differences of group scatters, a value of 1 specifies the strongest restriction.
#'
#' @return A list with entries:
#' \describe{
#'  \item{cluster}{A numerical vector of size n containing the cluster assignment for each observation. Cluster names are integer numbers from 1 to k, 0 indicates trimmed observations.}
#'  \item{n_clus}{Number of clusters actually found.}
#'  \item{obj}{he value of the objective function of the best (returned) solution.}
#' }
#'
#' @examples
#' x <- rbind(matrix(rnorm(100), ncol = 2), matrix(rnorm(100) + 2, ncol = 2),
#'         matrix(rnorm(100) + 4, ncol = 2))
#' ## robust cluster obtention from a sample x asking for 3 clusters,
#' ## trimming level 0.05 and constrain level 12
#' k <- 3; alpha <- 0.05; restr.fact <- 12
#' output = tclust_H(x = x, k = k, alpha = alpha, nstart = 50, iter.max = 20,
#'                  restr = 'eigen', restr.fact = restr.fact, sol_ini_p = FALSE, sol_ini = NA,
#'                  equal.weights = FALSE, trace = 0, zero.tol = 1e-16)
#' ## cluster assigment
#' output2 <- tclustWithInitialization(data.frame(x, output$cluster), x, 'points', 0.05, 10)
#'
#' @export
#'
tclustWithInitialization <- function(initialization, cytometry, i.sol.type = "points", trimming = 0.05, restr.fact = 1000) {
    t0 <- Sys.time()
    dim.cyto <- dim(cytometry)[2]
    if (i.sol.type == "points") {
        cosa <- data.frame(initialization[, 1:dim.cyto], cluster = dplyr::pull(initialization, (dim.cyto + 1)))
        n_clus <- length(levels(cosa$cluster))
        levels(cosa$cluster) <- 1:n_clus

        cosa_no_noise <- cosa[which(as.double(cosa$cluster) > 0), ]

        sol_inicial_0 <- lapply(1:n_clus, estimCovCellGeneral, cosa[, 1:dim.cyto], cosa$cluster)
        not_valid_estim <- which(is.na(sol_inicial_0))
        if (length(not_valid_estim) > 0) {
            sol_inicial_0 <- sol_inicial_0[-not_valid_estim]
        } else {
            sol_inicial_0 <- sol_inicial_0
        }
        weight_citos <- unlist(lapply(sol_inicial_0, function(x) x$weight))
        means_citos <- matrix(unlist(lapply(sol_inicial_0, function(x) x$mean)), ncol = length(weight_citos), byrow = FALSE)
        cov_citos <- array(dim = c(dim.cyto, dim.cyto, length(weight_citos)))
        for (i in 1:length(weight_citos)) {
            cov_citos[, , i] <- sol_inicial_0[[i]]$cov
        }

        sol_inicial <- list()
        sol_inicial$weights <- weight_citos
        sol_inicial$cov <- cov_citos
        sol_inicial$centers <- means_citos
        KK <- length(sol_inicial$weights)
        print(paste("tclust looking for k = ", KK, sep = ""))
        trimming <- 0.5 * exp(-(KK - 1)/5)
    } else {
        if (i.sol.type == "barycenter") {
            sol_inicial <- initialization
            KK <- length(sol_inicial$weights)
            print(paste("tclust looking for k = ", KK, sep = ""))
        }
    }

    t00 <- Sys.time()
    solution.tclust <- optimalFlow::tclust_H(as.matrix(cytometry[, 1:dim.cyto]), k = KK, alpha = trimming, restr.fact = restr.fact,
        nstart = 10, iter.max = 100, equal.weights = FALSE, sol_ini = sol_inicial, sol_ini_p = TRUE, restr = "eigen", trace = 0,
        zero.tol = 1e-16)
    t11 <- Sys.time()
    t11 - t00

    t1 <- Sys.time()
    n_clus <- length(table(solution.tclust$cluster)) - 1
    print(paste("tclust found k = ", n_clus, sep = ""))
    print(paste(difftime(t1, t0), units(difftime(t1, t0))))
    cat("\n")
    return(list(cluster = solution.tclust$cluster, n_clus = n_clus, obj = solution.tclust$obj))
}
