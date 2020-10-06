#' optimalFlowTemplates
#'
#' Returns a partition of the input clusterings with a respective consensus clustering for every group.
#'
#' @param database A list where each entry is a partition (clustering) represented as dataframe, of the same dimensions, where the last variable represents the labels of the partition.
#' @param database.names Names of the elements in the database.
#' @param cov.estimation How to estimate covariance matrices in each cluster of a partition. 'standard' is for using cov(), while 'robust' is for using robustbase::covMcd.
#' @param alpha.cov Only when cov.estimation = 'robust'. Indicates the value of alpha in robustbase::covMcd.
#' @param equal.weights.template If True, weights assigned to every cluster in a partion are uniform (1/number of clusters). If False, weights assigned to clusters are the proportions of points in every cluster compared to the total amount of points in the partition.
#' @param hclust.method Indicates what kind of hierarchical clustering to  do with the similarity distances matrix of the partitions. Takes values in c('complete', 'single', 'average', 'hdbscan', 'dbscan').
#' @param trimm.template Logical value. Indicates if it is allowed to not take into account some of the entries of database. Default is False.
#' @param templates.number Only if hclust.method in c('complete', 'single', 'average'). Indicates the number of clusters to use with cutree. If set to NA (default), plots the hierarchical tree and asks the user to introduce an appropriate number of clusters.
#' @param minPts Only if hclust.method in c('hdbscan', 'dbscan'). Indicates the value of argument minPts in dbscan::dbscan and dbscan::hdbscan.
#' @param eps Only if hclust.method = 'dbscan'. Indicates the value of eps in dbscan::dbscan.
#' @param consensus.method Sets the way of doing consensus clustering when clusters are viewed as Multivariate Distributions. Can take values in c('pooling', 'k-barycenter', 'hierarchical'). See details.
#' @param barycenters.number Only if consensus.method = 'k-barycenter'. Sets the number, k, of barycenters when using k-barycenters.
#' @param bar.repetitions Only if consensus.method = 'k-barycenter'. How many times to repeat the k-barycenters procedure. Equivalent to nstart in kmeans.
#' @param alpha.bar Only if consensus.method = 'k-barycenter'. The level of trimming allowed during the k-barycenters procedure.
#' @param bar.ini.method Only if consensus.method = 'k-barycenter'. Takes values in c('rnd', 'plus-plus'). See details.
#' @param consensus.minPts Only if consensus.method = 'hierarchical'. The value of argument minPts for dbscan::hdbscan.
#' @param cl.paral Number of cores to be used in parallel procedures.
#'
#' @return A list containting:
#' \describe{
#'  \item{templates}{A list representing the consensus clusterings for every group in the partition of the database. Each element of the list is a template partition. Hence it is a list itself, containig the cell types in the prototype, where each element has components: mean, cov, weight and type.}
#'  \item{clustering}{Clustering of the input partitions.}
#'  \item{database.elliptical}{A list containig each cytometry in the database viewed as a mixture distribution. Each element of the list is a cytometry viewed as a mixture. Hence it is a list itself, containig the cell types in the cytometry, where each element has components: mean, cov, weight and type.}
#' }
#'
#' @examples
#' # # We construct a simple database selecting only some of the Cytometries and some cell types for simplicity and for a better visualisation.
#' database <- buildDatabase(
#'   dataset_names = paste0('Cytometry', c(2:5, 7:9, 12:17, 19, 21)),
#'     population_ids = c('Monocytes', 'CD4+CD8-', 'Mature SIg Kappa', 'TCRgd-'))
#'
#' # # To select the appropriate number of templates, via hierarchical tree, in an interactive fashion and produce a clustering we can also use:
#' # templates.optimalFlow <- optimalFlowTemplates(database = database)
#'
#' templates.optimalFlow <- optimalFlowTemplates(database = database, templates.number = 5,
#'                                              cl.paral = 1)
#'
#' @references E del Barrio, H Inouzhe, JM Loubes, C Matran and A Mayo-Iscar. (2019) optimalFlow: Optimal-transport approach to flow cytometry gating and population matching. arXiv:1907.08006
#'
#' @import dplyr
#' @import optimalFlowData
#' @import rlang
#' @importFrom foreach %dopar%
#' @importFrom foreach %:%
#' @importFrom stats as.dist
#' @importFrom stats cutree
#' @importFrom stats cov
#' @importFrom stats hclust
#' @importFrom stats runif
#'
#' @export
#'
optimalFlowTemplates <- function(database, database.names = NULL, cov.estimation = "standard", alpha.cov = 0.85, equal.weights.template = TRUE,
    hclust.method = "complete", trimm.template = FALSE, templates.number = NA, minPts = 2, eps = 1, consensus.method = "pooling",
    barycenters.number = NA, bar.repetitions = 40, alpha.bar = 0.05, bar.ini.method = "plus-plus", consensus.minPts = 3, cl.paral = 1) {
    t0 <- Sys.time()
    op.syst <- .Platform$OS.type

    ####################### Step 0: Checking that entries are well defined #######################

    cov.estimation <- match.arg(cov.estimation, c("standard", "robust"))
    if (!is.double(alpha.cov)) {
        stop("alha.cov is not a doble")
    } else {
        if (alpha.cov < 0 | alpha.cov >= 1) {
            stop("alpha.cov is not in the range (0,1).")
        }
    }
    hclust.method <- match.arg(hclust.method, c("complete", "single", "average", "hdbscan", "dbscan"))
    if (!is.na(templates.number)) {
        if (!is.double(templates.number) & !is.integer(templates.number)) {
            stop("templates.number is not well defined.")
        } else {
            if (templates.number < 1) {
                stop("templates.number should be >= 1.")
            }
        }
    }
    if (hclust.method %in% c("hdbscan", "dbscan")) {
        minPts <- as.integer(minPts)
        if (is.na(minPts) | is.na(eps)) {
            stop("minPts or eps not well defined.")
        }
    }
    consensus.method <- match.arg(consensus.method, c("pooling", "k-barycenter", "hierarchical"))
    if (consensus.method == "k-barycenter") {
        barycenters.number <- as.integer(barycenters.number)
        bar.repetitions <- as.integer(bar.repetitions)
        if (is.na(bar.repetitions)) {
            stop("bar.repetitions must be an integer >= 1.")
        }
        if (!is.double(alpha.bar)) {
            stop("alha.bar is not a doble")
        } else {
            if (alpha.bar < 0 | alpha.bar >= 1) {
                stop("alpha.bar is not in the range [0,1).")
            }
        }
        bar.ini.method <- match.arg(bar.ini.method, c("rnd", "plus-plus"))
    } else {
        if (consensus.method == "hierarchical") {
            consensus.minPts <- as.integer(consensus.minPts)
            if (is.na(consensus.minPts)) {
                stop("consensus.minPts or consensus.eps not well defined.")
            }
        }
    }
    op.syst <- match.arg(op.syst, c("unix", "windows"))

    ####################### Step 1: Similarity distance based clustering #######################

    sys.time.step1.0 <- Sys.time()
    database.elliptical <- list()
    for (i in 1:length(database)) {
        dim.cyto <- dim(database[[i]])[2]
        database.elliptical[[i]] <- lapply(names(table(database[[i]][, dim.cyto])), estimCovCellGeneral, cytometry = database[[i]][,
            1:(dim.cyto - 1)], labels = database[[i]][, dim.cyto], type = cov.estimation, alpha = alpha.cov)
        database.elliptical[[i]] <- database.elliptical[[i]][!is.na(database.elliptical[[i]])]
    }

    # cl <- parallel::makePSOCKcluster(cl.paral) doParallel::registerDoParallel(cl) database.elliptical = foreach::foreach(i =
    # 1:length(database)) %dopar% { source('estimCovCellGeneral.R') database.elliptical0 = lapply(names(table(database[[i]][,
    # dim.cyto])), estimCovCellGeneral, cytometry = database[[i]][, 1:(dim.cyto - 1)], labels = database[[i]][, dim.cyto], type =
    # cov.estimation, alpha = alpha.cov) database.elliptical0 = database.elliptical0[!is.na(database.elliptical0)] }
    # parallel::stopCluster(cl)

    dim.cytos <- length(database.elliptical)
    if (!is.na(templates.number) & dim.cytos < templates.number){
      stop(paste("templates.number, currently ", templates.number,", should be betwen 1 and ", dim.cytos," (the number of elements in the database).", sep = ""))
    }
    wasser.cost <- array(0, dim = c(dim.cytos, dim.cytos))
    cytometries_elipses <- list()
    if ((dim.cytos - 2) * (dim.cytos - 2 + 1)/2 < cl.paral) {
        n.cores <- (dim.cytos - 2) * (dim.cytos - 2 + 1)/2
    } else {
        n.cores <- cl.paral
    }
    cl <- parallel::makePSOCKcluster(n.cores, setup_strategy = "sequential")
    doParallel::registerDoParallel(cl)
    transport_costs_list <- foreach::foreach(j = 1:(dim.cytos - 1)) %:% foreach::foreach(i = (j + 1):dim.cytos) %dopar% {
        optimalFlow::wasserCostFunction(j, i, database.elliptical, equal.weights.template)
    }
    wasser.cost[lower.tri(wasser.cost, diag = FALSE)] <- unlist(transport_costs_list)
    parallel::stopCluster(cl)


    if (hclust.method %in% c("complete", "single", "average")) {
        cytos.hcl <- hclust(as.dist(wasser.cost), method = hclust.method)
        if (is.na(templates.number)) {
            if (is.null(database.names)) {
                graphics::plot(cytos.hcl, main = NULL, ylab = NULL, xlab = "")
            } else {
                graphics::plot(cytos.hcl, labels = database.names, main = NULL, ylab = NULL, xlab = "")
            }
            templates.number <- as.integer(readline(prompt = "How many clusters should I look for : "))
            if (!is.na(templates.number) & dim.cytos < templates.number){
              stop(paste("templates.number, currently ", templates.number,", should be betwen 1 and ", dim.cytos," (the number of elements in the database).", sep = ""))
            }
            cytos.cluster <- cutree(cytos.hcl, k = templates.number)
        } else {
            cytos.cluster <- cutree(cytos.hcl, k = templates.number)
        }
    } else {
        if (hclust.method == "hdbscan") {
            cytos.cluster <- dbscan::hdbscan(as.dist(wasser.cost), minPts = minPts)$cluster
        } else {
            cytos.cluster <- dbscan::dbscan(as.dist(wasser.cost), minPts = minPts, eps = eps)$cluster
        }
    }
    if (!trimm.template) {
        indexes.trimmed <- which(cytos.cluster == 0)
        if (length(indexes.trimmed) > 0) {
            cytos.cluster <- cytos.cluster + length(indexes.trimmed)
            cytos.cluster[indexes.trimmed] <- seq(1, length(indexes.trimmed), 1)
        }
    }
    sys.time.step1.1 <- Sys.time()
    time.dif.1 <- difftime(sys.time.step1.1, sys.time.step1.0)
    print(paste("step 1:", time.dif.1, units(time.dif.1), sep = " "))

    ####################### Step 2: Consensus clustering #######################

    sys.time.step2.0 <- Sys.time()
    templates <- list()
    templates.cell.types <- list()
    if (consensus.method == "pooling") {
        ########## pooling clusters with the same labels and taking barycenter #########
        j <- 0
        for (i in sort(unique(cytos.cluster))) {
            if (i == 0) {
                next
            } else {
                j <- j + 1
                cell.types <- lapply(database.elliptical[cytos.cluster == i], function(x) lapply(x, function(y) y$type))
                cell.types <- sort(unique(unlist(cell.types)))
                templates.cell.types[[j]] <- cell.types
                if (length(cell.types) < cl.paral) {
                  n.cores <- length(cell.types)
                } else {
                  n.cores <- cl.paral
                }
                if (op.syst == "unix") {
                  templates[[j]] <- parallel::mclapply(cell.types, optimalFlow::estimationCellBarycenter, database.elliptical[cytos.cluster ==
                    i], mc.cores = n.cores)
                } else {
                  cl <- parallel::makePSOCKcluster(n.cores)
                  templates[[j]] <- parallel::parLapply(cl, cell.types, optimalFlow::estimationCellBarycenter, database.elliptical[cytos.cluster ==
                    i])
                  parallel::stopCluster(cl)
                }
            }
        }
    } else {
        ######## pooling all clusters for the partitions in the same group #######
        pooled.clustered.cytometries <- list()
        for (i in sort(unique(cytos.cluster))) {
            pooled.clustered.cytometries.0 <- list()
            l <- 0
            for (training.cytometry in database.elliptical[cytos.cluster == i]) {
                for (k in 1:length(training.cytometry)) {
                  l <- l + 1
                  pooled.clustered.cytometries.0[[l]] <- training.cytometry[[k]]
                }
            }
            pooled.clustered.cytometries[[as.character(i)]] <- pooled.clustered.cytometries.0
        }
        if (consensus.method == "hierarchical") {
            ######## Density based hierarchical clustering based on the wasserstein distance between clusters. For each resulting grouping we take
            ######## the one-barycenter.
            j <- 0
            for (i in sort(unique(cytos.cluster))) {
                if (i == 0) {
                  next
                } else {
                  j <- j + 1
                  if (sum(cytos.cluster == i) == 1) {
                    cell.types <- lapply(database.elliptical[cytos.cluster == i], function(x) lapply(x, function(y) y$type))
                    cell.types <- sort(unique(unlist(cell.types)))
                    template.partition <- cell.types
                  } else {
                    dim.cytos <- length(pooled.clustered.cytometries[[as.character(i)]])
                    wasser.cost.list <- wasser.cost <- array(0, dim = c(dim.cytos, dim.cytos))
                    if ((dim.cytos - 2) * (dim.cytos - 2 + 1)/2 < cl.paral) {
                      n.cores <- (dim.cytos - 2) * (dim.cytos - 2 + 1)/2
                    } else {
                      n.cores <- cl.paral
                    }
                    cl <- parallel::makePSOCKcluster(n.cores, setup_strategy = "sequential")
                    doParallel::registerDoParallel(cl)
                    wasser.cost.list <- foreach::foreach(k = 1:(dim.cytos - 1)) %:% foreach::foreach(l = (k + 1):dim.cytos) %dopar%
                      {
                        optimalFlow::w2dist(pooled.clustered.cytometries[[as.character(i)]][[k]], pooled.clustered.cytometries[[as.character(i)]][[l]])
                      }
                    wasser.cost[lower.tri(wasser.cost, diag = FALSE)] <- unlist(wasser.cost.list)
                    parallel::stopCluster(cl)
                    template.partition <- dbscan::hdbscan(as.dist(wasser.cost), minPts = consensus.minPts)$cluster
                  }
                  template.partition.labels <- sort(unique(template.partition))
                  if (length(which(template.partition.labels == 0)) > 0) {
                    template.partition.labels <- template.partition.labels[-1]
                  }
                  cl <- parallel::makePSOCKcluster(cl.paral, setup_strategy = "sequential")
                  doParallel::registerDoParallel(cl)
                  template <- foreach::foreach(k = template.partition.labels) %dopar% {
                    pooled.elements <- pooled.clustered.cytometries[[as.character(i)]][template.partition == k]
                    pooled.elements.type.proportion <- unlist(lapply(pooled.elements, function(x) x$type))
                    pooled.elements.type.proportion <- table(pooled.elements.type.proportion)
                    pooled.elements.type.proportion <- sort(pooled.elements.type.proportion/sum(pooled.elements.type.proportion),
                      decreasing = TRUE)
                    pooled.elements.covs <- lapply(pooled.elements, function(x) x$cov)
                    pooled.elements.weighted.means <- matrix(unlist(lapply(pooled.elements, function(x) x$mean)), ncol = length(pooled.elements),
                      byrow = FALSE)
                    pooled.elements.weights <- unlist(lapply(pooled.elements, function(x) x$weight))
                    mean.weight <- mean(pooled.elements.weights)
                    pooled.elements.weights <- pooled.elements.weights/sum(pooled.elements.weights)
                    pooled.mean <- colSums(t(pooled.elements.weighted.means) * pooled.elements.weights)
                    template.cell <- list(mean = pooled.mean, cov = optimalFlow:::GaussianBarycenters(pooled.elements.covs, pooled.elements.weights)$Barycenter,
                      weight = mean.weight, type = k, type.proportions = pooled.elements.type.proportion)
                  }
                  parallel::stopCluster(cl)
                  templates[[j]] <- template
                  templates.cell.types[[j]] <- template.partition.labels
                }
            }
        } else {
            j <- 0
            for (i in sort(unique(cytos.cluster))) {
                if (i == 0) {
                  next
                } else {
                  j <- j + 1
                  if (sum(cytos.cluster == i) == 1) {
                    cell.types <- lapply(database.elliptical[cytos.cluster == i], function(x) lapply(x, function(y) y$type))
                    cell.types <- sort(unique(unlist(cell.types)))
                    template.partition <- cell.types
                    template.partition.labels <- sort(unique(template.partition))
                    if (length(which(template.partition.labels == 0)) > 0) {
                      template.partition.labels <- template.partition.labels[-1]
                    }
                    cl <- parallel::makePSOCKcluster(cl.paral, setup_strategy = "sequential")
                    doParallel::registerDoParallel(cl)
                    template <- foreach::foreach(k = template.partition.labels) %dopar% {
                      pooled.elements <- pooled.clustered.cytometries[[as.character(i)]][template.partition == k]
                      pooled.elements.type.proportion <- unlist(lapply(pooled.elements, function(x) x$type))
                      pooled.elements.type.proportion <- table(pooled.elements.type.proportion)
                      pooled.elements.type.proportion <- sort(pooled.elements.type.proportion/sum(pooled.elements.type.proportion),
                        decreasing = TRUE)
                      pooled.elements.covs <- lapply(pooled.elements, function(x) x$cov)
                      pooled.elements.weighted.means <- matrix(unlist(lapply(pooled.elements, function(x) x$mean)), ncol = length(pooled.elements),
                        byrow = FALSE)
                      pooled.elements.weights <- unlist(lapply(pooled.elements, function(x) x$weight))
                      mean.weight <- mean(pooled.elements.weights)
                      pooled.elements.weights <- pooled.elements.weights/sum(pooled.elements.weights)
                      pooled.mean <- colSums(t(pooled.elements.weighted.means) * pooled.elements.weights)
                      template.cell <- list(mean = pooled.mean, cov = optimalFlow:::GaussianBarycenters(pooled.elements.covs, pooled.elements.weights)$Barycenter,
                        weight = mean.weight, type = k, type.proportions = pooled.elements.type.proportion)
                    }
                    parallel::stopCluster(cl)
                    templates[[j]] <- template
                    templates.cell.types[[j]] <- template.partition.labels
                  } else {
                    if (bar.repetitions < cl.paral) {
                      n.cores <- bar.repetitions
                    } else {
                      n.cores <- cl.paral
                    }
                    if (op.syst == "unix") {
                      if (is.na(barycenters.number)) {
                        barycenters.number.0 <- ceiling(mean(unlist(lapply(which(cytos.cluster == i), function(ii) length(database.elliptical[[ii]])))))
                        barycenters.number.0 <- round((1 + 0.1) * barycenters.number.0)
                        k.barycenter.list <- parallel::mclapply(1:bar.repetitions, optimalFlow:::wassersteinKBarycenter, k = barycenters.number.0,
                          alpha = alpha.bar, initialization = "plus-plus", pooled.clustered.cytometries[[as.character(i)]], mc.cores = n.cores)  #n.cores
                      } else {
                        if (barycenters.number > length(pooled.clustered.cytometries[[as.character(i)]])) {
                          stop(paste("barycenters.number, currently ", barycenters.number, ", should be between 1 and ", length(pooled.clustered.cytometries[[as.character(i)]]), " (the total number of elements in the group).", sep = ""))
                        }
                        k.barycenter.list <- parallel::mclapply(1:bar.repetitions, optimalFlow:::wassersteinKBarycenter, k = barycenters.number,
                          alpha = alpha.bar, initialization = "plus-plus", pooled.clustered.cytometries[[as.character(i)]], mc.cores = n.cores)  #n.cores
                      }
                    } else {
                      if (is.na(barycenters.number)) {
                        barycenters.number.0 <- ceiling(mean(unlist(lapply(which(cytos.cluster == i), function(ii) length(database.elliptical[[ii]])))))
                        barycenters.number.0 <- round((1 + 0.1) * barycenters.number.0)
                        cl <- parallel::makePSOCKcluster(n.cores, setup_strategy = "sequential")
                        k.barycenter.list <- parallel::parLapply(cl, 1:bar.repetitions, optimalFlow:::wassersteinKBarycenter, k = barycenters.number.0,
                          alpha = alpha.bar, initialization = "plus-plus", pooled.clustered.cytometries[[as.character(i)]])
                        parallel::stopCluster(cl)
                      } else {
                        if (barycenters.number > length(pooled.clustered.cytometries[[as.character(i)]])) {
                          stop(paste("barycenters.number, currently ", barycenters.number, ", should be between 1 and ", length(pooled.clustered.cytometries[[as.character(i)]]), " (the total number of elements in the group).", sep = ""))
                        }
                        cl <- parallel::makePSOCKcluster(n.cores, setup_strategy = "sequential")
                        k.barycenter.list <- parallel::parLapply(cl, 1:bar.repetitions, optimalFlow:::wassersteinKBarycenter, k = barycenters.number,
                          alpha = alpha.bar, initialization = "plus-plus", pooled.clustered.cytometries[[as.character(i)]])
                        parallel::stopCluster(cl)
                      }
                    }
                    wasser.var <- lapply(1:bar.repetitions, function(ii) k.barycenter.list[[ii]]$wasserstein.var)
                    if (min(unlist(wasser.var)) == Inf) {
                      stop("All executions of k-barycenter were unsuccesfull")
                    } else {
                      k.barycenter <- k.barycenter.list[[which.min(wasser.var)]]$wasserstein.k.barycenter
                      template.partition <- k.barycenter$cluster
                      template.partition.labels <- sort(unique(template.partition))
                      if (length(which(template.partition.labels == 0)) > 0) {
                        template.partition.labels <- template.partition.labels[-1]
                      }

                      if (length(template.partition.labels) < cl.paral) {
                        n.cores <- length(template.partition.labels)
                      } else {
                        n.cores <- cl.paral
                      }
                      if (op.syst == "unix") {
                        template <- parallel::mclapply(template.partition.labels, function(k) {
                          pooled.elements <- pooled.clustered.cytometries[[as.character(i)]][template.partition == k]
                          pooled.elements.type.proportion <- unlist(lapply(pooled.elements, function(x) x$type))
                          pooled.elements.type.proportion <- table(pooled.elements.type.proportion)
                          pooled.elements.type.proportion <- sort(pooled.elements.type.proportion/sum(pooled.elements.type.proportion),
                            decreasing = TRUE)
                          pooled.elements.weights <- unlist(lapply(pooled.elements, function(x) x$weight))
                          mean.weight <- mean(pooled.elements.weights)
                          template.cell <- list(mean = k.barycenter$baricentro[[which(template.partition.labels == k)]]$mean, cov = k.barycenter$baricentro[[which(template.partition.labels ==
                            k)]]$cov, weight = mean.weight, type = k, type.proportions = pooled.elements.type.proportion)
                        }, mc.cores = n.cores)
                      } else {
                        cl <- parallel::makePSOCKcluster(n.cores)
                        template <- parallel::parLapply(cl, template.partition.labels, function(k) {
                          pooled.elements <- pooled.clustered.cytometries[[as.character(i)]][template.partition == k]
                          pooled.elements.type.proportion <- unlist(lapply(pooled.elements, function(x) x$type))
                          pooled.elements.type.proportion <- table(pooled.elements.type.proportion)
                          pooled.elements.type.proportion <- sort(pooled.elements.type.proportion/sum(pooled.elements.type.proportion),
                            decreasing = TRUE)
                          pooled.elements.weights <- unlist(lapply(pooled.elements, function(x) x$weight))
                          mean.weight <- mean(pooled.elements.weights)
                          template.cell <- list(mean = k.barycenter$baricentro[[which(template.partition.labels == k)]]$mean, cov = k.barycenter$baricentro[[which(template.partition.labels ==
                            k)]]$cov, weight = mean.weight, type = k, type.proportions = pooled.elements.type.proportion)
                        })
                        parallel::stopCluster(cl)
                      }
                      templates[[j]] <- template
                      templates.cell.types[[j]] <- template.partition.labels
                    }
                  }
                }
            }
        }
    }
    sys.time.step2.1 <- Sys.time()
    time.dif.2 <- difftime(sys.time.step2.1, sys.time.step2.0)
    print(paste("step 2:", time.dif.2, units(time.dif.2), sep = " "))
    t1 <- Sys.time()
    time.dif.total <- difftime(sys.time.step2.1, t0)
    print(paste("Execution time:", time.dif.total, units(time.dif.total), sep = " "))
    return(list(templates = templates, clustering = cytos.cluster, database.elliptical = database.elliptical))
}
