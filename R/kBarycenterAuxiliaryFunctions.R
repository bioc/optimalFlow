#' distGaussianMean
#'
#' Computes the componente relative to the means of the squared 2-Wasserstein distance.
#'
#' @param N1 A multivariate normal distribution as list with elements mean and cov.
#' @param N2 A multivariate normal distribution as list with elements mean and cov.
#'
#' @return A double equivalent to the squared Euclidean distance between the means.
#'
#' @examples
#' distGaussianMean(list(mean = c(-1, -1), cov = diag(1, 2)), list(mean = c(1, 1), cov = diag(1, 2)))
#'
#' @noRd
#'
distGaussianMean <- function(N1, N2) {
    di.mean <- sum((N1$mean - N2$mean)^2)
    distancia <- di.mean
    distancia
}
#' distGaussianCov
#'
#' Computes the componente relative to the covariances of the squared 2-Wasserstein distance.
#'
#' @param N1 A multivariate normal distribution as list with elements mean and cov.
#' @param N2 A multivariate normal distribution as list with elements mean and cov.
#'
#' @return A double representing the squared 2-Wassersteindistance between normals with the same mean.
#'
#' @examples
#' distGaussianCov(list(mean = c(1, 1), cov = diag(2, 2)), list(mean = c(1, 1), cov = diag(1, 2)))
#'
#' @noRd
#'
distGaussianCov <- function(N1, N2) {
    A <- N1$cov
    B <- N2$cov
    sqrt.matrix <- function(C) {
        e <- eigen(C)
        e$values <- abs(e$values)
        sqrtA <- e$vectors %*% diag(sqrt(e$values)) %*% t(e$vectors)
        sqrtA
    }
    sqrtA <- sqrt.matrix(A)
    M1 <- sqrtA %*% B %*% sqrtA
    M2 <- sqrt.matrix(M1)
    M3 <- A + B - 2 * M2
    di.cov <- sum(diag(M3))
    distancia <- di.cov
    distancia
}
#' distGaussian
#'
#' Computes the squared 2-Wasserstein distance.
#'
#' @param N1 A multivariate normal distribution as list with elements mean and cov.
#' @param N2 A multivariate normal distribution as list with elements mean and cov.
#'
#' @return A double representing the squared 2-Wasserstein distance between normals.
#'
#' @examples
#' distGaussian(list(mean = c(-1, -1), cov = diag(2, 2)), list(mean = c(1, 1), cov = diag(1, 2)))
#'
#' @noRd
#'
distGaussian <- function(N1, N2) {
    di.mean <- sum((N1$mean - N2$mean)^2)
    A <- N1$cov
    B <- N2$cov
    sqrt.matrix <- function(C) {
        e <- eigen(C)
        e$values <- abs(e$values)
        sqrtA <- e$vectors %*% diag(sqrt(e$values)) %*% t(e$vectors)
        sqrtA
    }
    sqrtA <- sqrt.matrix(A)
    M1 <- sqrtA %*% B %*% sqrtA
    M2 <- sqrt.matrix(M1)
    M3 <- A + B - 2 * M2
    di.cov <- sum(diag(M3))
    distancia <- di.mean + di.cov
    distancia
}
#' GaussianBarycenters
#'
#' Calculates the (1-)barrycenter of a mixture of normals with the same mean.
#'
#' @param matrices A list of covariance matrices.
#' @param weight A vector of weights associated to the covariance matrices.
#'
#' @return A list formed by:
#' \describe{
#'  \item{Barycenter}{Returns the barycenter, i.e., a covariance matrix.}
#'  \item{Variation}{A double representing the Wasserstein Variation.}
#'  \item{Num.iter}{The number of iterations to achieve the stopping criteria.}
#' }
#'
#' @examples
#' GaussianBarycenters(list(diag(2, 2),diag(1, 2)), c(0.5, 0.5))
#'
#' @noRd
#'
GaussianBarycenters <- function(matrices, weight) {
    d <- ncol(matrices[[1]])
    k <- length(weight)
    Sn <- diag(d)
    trfix <- 0
    for (j in 1:k) trfix <- trfix + sum(weight[j] * sum(diag(matrices[[j]])))
    main.iter <- function(S) {
        VS <- sum(diag(S)) + trfix
        aux1 <- eigen(S)
        aux1$values <- abs(aux1$values)
        S12 <- aux1$vectors %*% diag(aux1$values^{
            1/2
        }) %*% t(aux1$vectors)
        S12inv <- aux1$vectors %*% diag(aux1$values^{
            -1/2
        }) %*% t(aux1$vectors)
        Saux <- list()
        trmen <- 0
        for (j in 1:k) {
            temp <- eigen(S12 %*% matrices[[j]] %*% S12)
            temp$values <- abs(temp$values)
            Saux[[j]] <- temp$vectors %*% diag(temp$values^{
                1/2
            }) %*% t(temp$vectors)
            Saux[[j]] <- weight[j] * Saux[[j]]
            trj <- sum(diag(Saux[[j]]))
            trmen <- trmen + trj
        }
        VS <- VS - 2 * trmen
        S <- Reduce("+", Saux)
        list(iterant = S, V = VS)
    }
    nuevo <- main.iter(Sn)
    Vold <- nuevo$V
    Sn <- nuevo$iterant
    tol <- 5
    n.iter <- 1
    while (tol > 1e-09) {
        nuevo <- main.iter(Sn)
        n.iter <- n.iter + 1
        Vnew <- nuevo$V
        tol <- Vold - Vnew
        Vold <- Vnew
        Sn <- nuevo$iterant
    }
    list(Barycenter = Sn, Variation = Vnew, Num.iter = n.iter)
}
#' kcenter
#'
#' Calculates the k-barycenter of a list of multivaraite normals for a given
#'  number of clusters ,k, and an assignation of each normal to the respective group
#'
#' @param points A list of multivariate normals, where each element is a list with values mean and cov.
#' @param kk The number k of groups for the k-barycenter.
#' @param center.asigned A vector indicating to which group each normal should be assigned.
#'
#' @return A list with values:
#' \describe{
#'  \item{kcenters}{A list of the elements of the k-barycenter, i.e., a list of k lists containing means and covariances.}
#'  \item{t.variation}{A double giving the trimmed wasserstein variation.}
#' }
#'
#' @examples
#' normals <- list(list(mean = c(1, 1), cov = diag(2, 2)), list(mean = c(1, 1), cov = diag(1, 2)),
#'  list(mean = c(3, 3), cov = diag(1, 2)))
#' kcenter(normals, 2, c(1, 1, 2))
#'
#' @noRd
#'
kcenter <- function(points, kk, center.asigned) {
    new.bar <- list()
    Trimmed.Variation.Prov <- 0
    for (a in 1:kk) {
        cluster.a <- points[which(center.asigned == a)]
        n.a <- length(cluster.a)
        means.a <- NULL
        covs.a <- list()
        for (i in 1:n.a) {
            means.a <- cbind(means.a, cluster.a[[i]]$mean)
            covs.a[[i]] <- cluster.a[[i]]$cov
        }
        mean.cluster.a <- apply(means.a, 1, mean)
        var.means.cluster.a <- mean(apply(means.a^2, 2, sum)) - sum(mean.cluster.a^2)
        Alg.Baricenter <- GaussianBarycenters(covs.a, weight = rep(1, n.a)/n.a)
        cov.cluster.a <- Alg.Baricenter$Barycenter
        new.bar[[a]] <- list(mean = mean.cluster.a, cov = cov.cluster.a)
        Trimmed.Variation.Prov <- Trimmed.Variation.Prov + Alg.Baricenter$Variation + var.means.cluster.a
    }
    resultado <- list(kcenters = new.bar, t.variation = Trimmed.Variation.Prov)
    resultado
}
#' wasserMinDist
#'
#' For two lists of multivariate normals calcualtes the closest member,
#'  in 2-Wasserstein distance, of the later list to each element of the former.
#'
#' @param points List of multivariate normals, where each element is a list with values mean and cov.
#' @param centres List of multivariate normals, where each element is a list with values mean and cov.
#'
#' @return A matrix where each column indicates the distance between the respective entry in points and
#'  the closest element in centers and the index of this closest element.
#'
#' @examples
#' normals <- list(list(mean = c(1, 1), cov = diag(2, 2)), list(mean = c(1, 1), cov = diag(1, 2)),
#'  list(mean = c(3, 3), cov = diag(1, 2)))
#' k_barycenter <- kcenter(normals, 2, c(1, 1, 2))$kcenters
#' wasserMinDist(normals, k_barycenter)
#'
#' @noRd
#'
wasserMinDist <- function(points, centres) {
    uu <- length(points)
    kk <- length(centres)
    aux.f1 <- function(u) {
        N1 <- points[[u]]
        aux.f2 <- function(j) distGaussian(N1, centres[[j]])
        distancias <- sapply(1:kk, aux.f2)
        resultado <- c(min(distancias), which.min(distancias))
        resultado
    }
    resultado <- sapply(1:uu, aux.f1)
    resultado
}
#' trimmedMinDist
#'
#' For two lists of multivaraite normals, points and centers, returns the index of which element
#'  in centers is closest to each element in points when some trimming is allowed.
#'
#' @param points List of multivariate normals, where each element is a list with values mean and cov.
#' @param centres List of multivariate normals, where each element is a list with values mean and cov.
#' @param alpha Level of triming.
#'
#' @return A vector with the index of which element in centers is closest to each element in points.
#'
#' @examples
#' normals <- list(list(mean = c(1, 1), cov = diag(2, 2)),list(mean = c(1, 1), cov = diag(1, 2)),
#'  list(mean = c(3, 3), cov = diag(1, 2)))
#' k_barycenter <- kcenter(normals, 2, c(1, 1, 2))$kcenters
#' trimmedMinDist(normals, k_barycenter, 0)
#'
#' @noRd
#'
trimmedMinDist <- function(points, centres, alpha = 0.1) {
    asignacion <- wasserMinDist(points, centres)
    orden.dista <- order(asignacion[1, ])
    threshold <- sort(asignacion[1, ])[floor(length(points) * (1 - alpha))]
    asignacion[2, (asignacion[1, ] > threshold)] <- 0
    resultado <- asignacion[2, ]
    resultado
}

