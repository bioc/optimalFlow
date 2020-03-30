#' trimmedKBarycenter
#'
#' Calculates a 2-Wasserstein k-barycenter of a list of multivariate normal distributions.
#'
#' @param k Number k of elements in the k-barycenter.
#' @param alpha0 Level of trimming.
#' @param type.ini of initialization in c('rnd', 'plus-plus'). 'rnd' makes the common random
#' initilaization while 'plus-plus' initializes in a similar fashion to k-means++.
#' @param reps.list List of multivariate normals for which the trimmed k-barycenter should be performed.
#'
#' @return A list with values:
#' \describe{
#'  \item{variacion_wasser}{A double giving the Waserstein variation.}
#'  \item{baricentro}{A list of k elements, each of which is a member of the k-barycenter.
#'   Each eement is a normal distribution characterized by a mean and a covariance.}
#'  \item{cluster}{The assignment of the original entries to each member of the k-barycenter.}
#' }
#'
#' @examples
#' normals <- list(list(mean = c(1, 1), cov = diag(2, 2)), list(mean = c(1, 1),cov = diag(1, 2)),
#'  list(mean = c(3, 3), cov = diag(1, 2)))
#' trimmedKBarycenter(2, 0, 'rnd', normals)
#'
#' @export
#'
trimmedKBarycenter <- function(k, alpha0, type.ini = "rnd", reps.list) {
    if (type.ini == "rnd") {
        indices_prov <- sample(1:(length(reps.list)), k, replace = FALSE)
        indices_prov <- sort(indices_prov)
        bar.prov <- reps.list[indices_prov]
    } else {
        if (type.ini == "plus-plus") {
            length_rep_lista <- length(reps.list)
            j <- 1
            c0 <- sample(1:(length_rep_lista), 1)
            centers <- list(reps.list[[c0]])
            
            if (k > 1) {
                cent_pos <- (1:(length_rep_lista))[-c0]
                dist_w2 <- array(0, dim = c(k, length_rep_lista))
                dist_w2[1, ] <- abs(unlist(lapply(reps.list, distGaussian, centers[[1]])))
                D2 <- dist_w2[1, ][-c0]/sum(dist_w2[1, ][-c0])
                
                for (j in 2:k) {
                  c0 <- sample(cent_pos, 1, prob = D2)
                  c0_pos <- which(cent_pos == c0)
                  cent_pos <- cent_pos[-c0_pos]
                  centers[[j]] <- reps.list[[c0]]
                  dist_w2[j, ] <- unlist(lapply(reps.list, distGaussian, centers[[j]]))
                  s <- 0
                  D2 <- array(0, dim = length(cent_pos))
                  for (i in cent_pos) {
                    s <- s + 1
                    d_min <- min(abs(dist_w2[1:j, i]))
                    D2[s] <- d_min
                  }
                  D2 <- D2/sum(D2)
                }
            }
            bar.prov <- centers
        } else {
            return("Inicialization not recognised")
        }
    }
    
    c.asig <- trimmedMinDist(reps.list, bar.prov, alpha = alpha0)
    update.step <- kcenter(reps.list, k, c.asig)
    new.bar <- update.step$kcenters
    Trimmed.Variation.Prov <- update.step$t.variation
    Trimmed.Variation.Old <- Trimmed.Variation.Prov
    bar.prov <- new.bar
    
    diferencia <- 10
    iter <- 0
    while (diferencia > 1e-07) {
        iter <- iter + 1
        c.asig <- trimmedMinDist(reps.list, bar.prov, alpha = alpha0)
        update.step <- kcenter(reps.list, k, c.asig)
        new.bar <- update.step$kcenters
        Trimmed.Variation.Prov <- update.step$t.variation
        diferencia <- Trimmed.Variation.Old - Trimmed.Variation.Prov
        Trimmed.Variation.Old <- Trimmed.Variation.Prov
        bar.prov <- new.bar
    }
    Trimmed.Variation <- Trimmed.Variation.Prov
    list(variacion_wasser = Trimmed.Variation, baricentro = bar.prov, cluster = c.asig)
}
