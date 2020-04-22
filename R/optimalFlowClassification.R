#' optimalFlowClassification
#'
#' Performs a supervised classification of input data when a database and a partition of the database are provided.
#'
#' @param X Datasample to be classified.
#' @param database A list where each entry is a partition (clustering) represented as dataframe, of the same dimensions, where the last variable represents the labels of the partition.
#' @param templates List of the consensus clusterings for every group in the partition of the database obtained by optimalFlowTemplates
#' @param consensus.method The consensus.method value that was used in optimalFlowTemplates.
#' @param cov.estimation How to estimate covariance matrices in each cluster of a partition. "standard" is for using cov(), while "robust" is for using robustbase::covMcd.
#' @param alpha.cov Only when cov.estimation = "robust". Indicates the value of alpha in robustbase::covMcd.
#' @param initial.method Indicates how to obtain a partition of X. Takes values in  c("supervized", "unsupervized"). Supervized uses tclust initilized by templates. Unsupevized usese flowMeans.
#' @param max.clusters The maximum numbers of clusters for flowMeans. Only when initial.method = unsupervized.
#' @param alpha.tclust Level of trimming allowed fo tclust. Only when initial.method = supervized.
#' @param restr.factor.tclust Fixes the restr.fact parameter in tclust. Only when initial.method = supervized.
#' @param classif.method Indicates what type of supervised learning we want to do. Takes values on c("matching", "qda", "random forest").
#' @param qda.bar Only if classif.method = "qda". If True then the appropriate consensus clustering (template, prototype) is used for learning. If False, the closest partition in the appropriate group is used.
#' @param cost.function Only if classif.method = "matching". Indicates the cost function, distance between clusters, to be used for label matching.
#' @param cl.paral Number of cores to be used in parallel procedures.
#' @param equal.weights.voting only when classif.method = "qda" and qda.bar =F, or when  classif.method = "random forest". Indicates the weights structure when looking for the most similar partition in a group.
#' @param equal.weights.template If True, weights assigned to every cluster in a partion are uniform (1/number of clusters). If False, weights assigned to clusters are the proportions of points in every cluster compared to the total amount of points in the partition.
#'
#' @return A list formed by:
#' \describe{
#' \item{cluster}{Labels assigned to the input data.}
#' \item{clusterings}{A list that contains the initial unsupervized or semi-supervized clusterings of the cytometry of interest. Can have as much entries as the number of templates in the semi-supervized case (initial.method = "supervized), or only one entry in the case of initial.method = "unsupervized". Each entry is a list where the most relevant argument for the clusterings is cluster.}
#' \item{assigned.template.index}{Label of the group for which the template is closer to the data. When classical qda or random forest ares used for classification there is a secon argument indicating the index of the cytometry in the cluster used for learning.}
#' \item{cluster.vote}{Only when classif.method = "matching" or when consensus.method in c("hierarchical", "k-barycenter"). Vote on the type of every label in the partition of the data. In essence, cluster + cluster.vote return a fuzzy clustering of the data of interest.}
#' }
#'
#' @examples
#' # # We construct a simple database selecting only some of the Cytometries and some cell types for simplicity and for a better visualisation.
#' database <- buildDatabase(
#'   dataset_names = paste0('Cytometry', c(2:5, 7:9, 12:17, 19, 21)),
#'     population_ids = c('Monocytes', 'CD4+CD8-', 'Mature SIg Kappa', 'TCRgd-'))
#' # # To select the appropriate number of templates, via hierarchical tree, in an interactive fashion and produce a clustering we can also use:
#' # templates.optimalFlow <- optimalFlowTemplates(database = database)
#' templates.optimalFlow <- optimalFlowTemplates(database = database, templates.number = 5,
#'                                              cl.paral = 1)
#' classification.optimalFlow <- optimalFlowClassification(Cytometry1[
#'   which(match(Cytometry1$`Population ID (name)`,c("Monocytes", "CD4+CD8-", "Mature SIg Kappa",
#'                                                   "TCRgd-"), nomatch = 0) > 0), 1:10], database, templates.optimalFlow, cl.paral = 1)
#' scoreF1.optimalFlow <- optimalFlow::f1Score(classification.optimalFlow$cluster,
#'                                            Cytometry1[which(match(Cytometry1$`Population ID (name)`,
#'                                                                                  c("Monocytes", "CD4+CD8-", "Mature SIg Kappa", "TCRgd-"), nomatch = 0) > 0),], noise.types)
#'
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
#' @importFrom stats predict
#' @import Rfast
#' @import flowMeans
#' @export
#'
optimalFlowClassification <- function (X, database, templates, consensus.method = "pooling", cov.estimation = "standard",
                                      alpha.cov = 0.85, initial.method = "supervized", max.clusters = NA, alpha.tclust = 0,
                                      restr.factor.tclust = 1000, classif.method = "qda", qda.bar = TRUE,
                                      cost.function = "points", cl.paral = 1, equal.weights.voting = TRUE,
                                      equal.weights.template = TRUE) {
    t0 <- Sys.time()
    op.syst <- .Platform$OS.type

    ####################### Step 0: Checking that entries are well defined #######################

    if (!is.matrix(X) & !is.data.frame(X)) {
        stop("Data is not a matrix or a dataframe.")
    }
    cov.estimation <- match.arg(cov.estimation, c("standard", "robust"))
    if (!is.double(alpha.cov)) {
        stop("alha.cov is not a doble")
    } else {
        if (alpha.cov < 0 | alpha.cov >= 1) {
            stop("alpha.cov is not in the range (0,1).")
        }
    }
    initial.method <- match.arg(initial.method, c("supervized", "unsupervized"))
    if (!is.double(alpha.tclust)) {
        stop("alha.tclust is not a doble")
    } else {
        if (alpha.tclust < 0 | alpha.tclust >= 1) {
            stop("alpha.tclust is not in the range [0,1).")
        }
    }
    classif.method <- match.arg(classif.method, c("matching", "qda", "random forest"))
    cost.function <- match.arg(cost.function, c("points", "ellipses"))
    op.syst <- match.arg(op.syst, c("unix", "windows"))

    ####################### Step 1: tclust step #######################

    sys.time.step1.0 <- Sys.time()
    dim.cyto <- dim(X)[2]
    if (initial.method == "unsupervized") {
        # stop('unsupervized initialization is not yet implemented.')
        assigned.tclust.index <- 1
        tclust.results <- list()
        tclust.results[[assigned.tclust.index]] <- list(cluster = flowMeans::flowMeans(X, MaxN = max.clusters, nstart = 10, Standardize = FALSE)@Label)

    } else {
        if (length(templates$templates) < cl.paral) {
            n.cores <- length(templates$templates)
        } else {
            n.cores <- cl.paral
        }
        if (op.syst == "unix") {
            tclust.results <- parallel::mclapply(templates$templates, function(template, alpha.tclust, restr.factor.tclust) {
                initial.weights <- unlist(lapply(template, function(x) x$weight))
                initial.means <- matrix(unlist(lapply(template, function(x) x$mean)), ncol = length(initial.weights), byrow = FALSE)

                initial.covs <- array(dim = c(dim.cyto, dim.cyto, length(initial.weights)))
                for (i in 1:length(initial.weights)) {
                  initial.covs[, , i] <- template[[i]]$cov
                }

                template.tclust <- list()
                template.tclust$cov <- initial.covs
                template.tclust$centers <- initial.means
                template.tclust$weights <- initial.weights/sum(initial.weights)
                clustering.tclust <- optimalFlow::tclustWithInitialization(template.tclust, X, i.sol.type = "barycenter", trimming = alpha.tclust,
                  restr.fact = restr.factor.tclust)
                return(clustering.tclust)
            }, alpha.tclust = alpha.tclust, restr.factor.tclust = restr.factor.tclust, mc.cores = n.cores)
        } else {
            cl <- parallel::makePSOCKcluster(n.cores)
            tclust.results <- parallel::parLapply(cl, templates$templates, function(template, alpha.tclust, restr.factor.tclust) {
                initial.weights <- unlist(lapply(template, function(x) x$weight))
                initial.means <- matrix(unlist(lapply(template, function(x) x$mean)), ncol = length(initial.weights), byrow = FALSE)

                initial.covs <- array(dim = c(dim.cyto, dim.cyto, length(initial.weights)))
                for (i in 1:length(initial.weights)) {
                  initial.covs[, , i] <- template[[i]]$cov
                }

                template.tclust <- list()
                template.tclust$cov <- initial.covs
                template.tclust$centers <- initial.means
                template.tclust$weights <- initial.weights/sum(initial.weights)
                clustering.tclust <- optimalFlow::tclustWithInitialization(template.tclust, X, i.sol.type = "barycenter", trimming = alpha.tclust,
                  restr.fact = restr.factor.tclust)
                return(clustering.tclust)
            }, alpha.tclust = alpha.tclust, restr.factor.tclust = restr.factor.tclust)
            parallel::stopCluster(cl)
        }
        objective.tclust <- unlist(lapply(tclust.results, function(x) x$obj))
        assigned.tclust.index <- which.max(objective.tclust)

    }

    sys.time.step1.1 <- Sys.time()
    time.dif.1 <- difftime(sys.time.step1.1, sys.time.step1.0)
    print(paste("step 1:", time.dif.1, units(time.dif.1), sep = " "))

    ####################### Step 2: Similarity distance assignation #######################

    sys.time.step2.0 <- Sys.time()

    data.elliptical <- lapply(sort(unique(tclust.results[[assigned.tclust.index]]$cluster)), optimalFlow::estimCovCellGeneral, cytometry = X,
        labels = tclust.results[[assigned.tclust.index]]$cluster, type = cov.estimation, alpha = alpha.cov)
    data.elliptical <- data.elliptical[!is.na(data.elliptical)]

    similarity.distances <- optimalFlow::costWasserMatchingEllipse(data.elliptical, templates$templates, equal.weights.template)
    print("Similarity distances to templates:")
    print(similarity.distances)
    assigned.template.index <- which.min(similarity.distances)
    assigned.template <- templates$templates[[assigned.template.index]]

    sys.time.step2.1 <- Sys.time()
    time.dif.2 <- difftime(sys.time.step2.1, sys.time.step2.0)
    print(paste("step 2:", time.dif.2, units(time.dif.2), sep = " "))

    ####################### Step 3: Labelling #######################

    sys.time.step3.0 <- Sys.time()

    if (classif.method == "qda" & qda.bar) {
        weight.groups <- unlist(lapply(assigned.template, function(x) x$weight))
        if (abs(sum(weight.groups) - 1) > 10^(-8)) {
            for (i in 1:length(assigned.template)) {
                assigned.template[[i]]$weight <- assigned.template[[i]]$weight/sum(weight.groups)
            }
        }
        if (length(assigned.template) < cl.paral) {
            n.cores <- length(assigned.template)
        } else {
            n.cores <- cl.paral
        }
        if (op.syst == "unix") {
            qda.assignation <- parallel::mclapply(assigned.template, optimalFlow::qdaClassification, data = X, mc.cores = n.cores)
        } else {
            cl <- parallel::makePSOCKcluster(n.cores)
            qda.assignation <- parallel::parLapply(cl, assigned.template, optimalFlow::qdaClassification, data = X)
            parallel::stopCluster(cl)
        }
        qda.assignation <- qda.assignation[!is.na(qda.assignation)]
        qda.assignation <- matrix(unlist(qda.assignation), ncol = length(qda.assignation), byrow = FALSE)
        qda.groups <- Rfast::rowMaxs(qda.assignation)
        unique.qda.groups <- sort(unique(qda.groups))
        qda.groups <- factor(qda.groups)
        levels(qda.groups) <- unlist(lapply(assigned.template, function(x) x$type))[!is.na(qda.assignation)]

        if (consensus.method == "hierarchical" | consensus.method == "k-barycenter") {
            vote <- list()
            if (sum(is.na(as.integer(levels(qda.groups)))) > 0) {
                range.qda.groups <- levels(qda.groups)
            } else {
                range.qda.groups <- as.integer(levels(qda.groups))
            }
            for (i in range.qda.groups) {
                vote[[as.character(i)]] <- data.frame(cell = names(assigned.template[[which(range.qda.groups == i)]]$type.proportions),
                  simple.proportion = as.vector(assigned.template[[which(range.qda.groups == i)]]$type.proportions))
            }
            sys.time.step3.1 <- Sys.time()
            time.dif.3 <- difftime(sys.time.step3.1, sys.time.step3.0)
            print(paste("step 3:", time.dif.3, units(time.dif.3), sep = " "))
            t1 <- Sys.time()
            print(difftime(t1, t0))
            return(list(cluster = qda.groups, cluster.vote = vote, vote.original.proportion = lapply(assigned.template, function(x) x$type.proportions),
                clusterings = tclust.results, assigned.template.index = assigned.template.index))
        } else {
            sys.time.step3.1 <- Sys.time()
            time.dif.3 <- difftime(sys.time.step3.1, sys.time.step3.0)
            print(paste("step 3:", time.dif.3, units(time.dif.3), sep = " "))
            t1 <- Sys.time()
            print(difftime(t1, t0))
            return(list(cluster = qda.groups, clusterings = tclust.results, assigned.template.index = assigned.template.index))
        }
    } else {
        if (classif.method == "matching") {
            if (cost.function == "points") {
                cytos.cluster <- templates$clustering
                if (sum(cytos.cluster == assigned.template.index) < cl.paral) {
                  n.cores <- sum(cytos.cluster == assigned.template.index)
                } else {
                  n.cores <- cl.paral
                }
                vote <- optimalFlow::voteLabelTransfer(type = cost.function, test.partition = tclust.results[[assigned.tclust.index]]$cluster,
                  test.cytometry = X, training.cytometries = database[cytos.cluster == assigned.template.index], test = 1, op.syst = op.syst,
                  cl.paral = n.cores, equal.weights = equal.weights.voting)
                sys.time.step3.1 <- Sys.time()
                time.dif.3 <- difftime(sys.time.step3.1, sys.time.step3.0)
                print(paste("step 3:", time.dif.3, units(time.dif.3), sep = " "))
                t1 <- Sys.time()
                print(difftime(t1, t0))
                return(list(cluster = tclust.results[[assigned.tclust.index]]$cluster, cluster.vote = vote$final.vote[[1]], clusterings = tclust.results,
                  assigned.template.index = assigned.template.index))

            } else {
                cytos.cluster <- templates$clustering
                vote <- optimalFlow::voteLabelTransfer(type = cost.function, test.partition.ellipse = data.elliptical, training.cytometries.barycenter = assigned.template,
                  test = 1, op.syst = op.syst, cl.paral = 1, equal.weights = equal.weights.voting)
                sys.time.step3.1 <- Sys.time()
                time.dif.3 <- difftime(sys.time.step3.1, sys.time.step3.0)
                print(paste("step 3:", time.dif.3, units(time.dif.3), sep = " "))
                t1 <- Sys.time()
                print(difftime(t1, t0))
                if (consensus.method == "hierarchical" | consensus.method == "k-barycenter") {
                  vote.original.proportion <- lapply(assigned.template, function(x) x$type.proportions)
                  if (sum(cytos.cluster == assigned.template.index) == 1) {
                    # print(vote) print(vote$final.vote[[1]]) final.vote = vote$fianl.vote[[1]] print(final.vote)
                    return(list(cluster = tclust.results[[assigned.tclust.index]]$cluster, cluster.vote = vote$final.vote[[1]],
                      clusterings = tclust.results, assigned.template.index = assigned.template.index))
                  } else {
                    final.vote <- voteTransformation(vote$final.vote[[1]], vote.original.proportion)
                    return(list(cluster = tclust.results[[assigned.tclust.index]]$cluster, cluster.vote = final.vote, clusterings = tclust.results,
                      assigned.template.index = assigned.template.index))
                  }
                } else {
                  return(list(cluster = tclust.results[[assigned.tclust.index]]$cluster, cluster.vote = vote$final.vote[[1]], clusterings = tclust.results,
                    assigned.template.index = assigned.template.index))
                }
            }
        } else {
            cytos.cluster <- templates$clustering
            database.elliptical <- templates$database.elliptical
            similarity.distances <- optimalFlow::costWasserMatchingEllipse(data.elliptical, database.elliptical[cytos.cluster ==
                assigned.template.index], equal.weights.template)
            assigned.training.cytometry.index <- which.min(similarity.distances)
            if (classif.method == "qda") {
                assigned.training.cytometry <- database.elliptical[cytos.cluster == assigned.template.index][[assigned.training.cytometry.index]]
                if (length(assigned.training.cytometry) < cl.paral) {
                  n.cores <- length(assigned.training.cytometry)
                } else {
                  n.cores <- cl.paral
                }
                if (op.syst == "unix") {
                  qda.assignation <- parallel::mclapply(assigned.training.cytometry, optimalFlow::qdaClassification, data = X, mc.cores = n.cores)
                } else {
                  cl <- parallel::makePSOCKcluster(n.cores)
                  qda.assignation <- parallel::parLapply(cl, assigned.training.cytometry, optimalFlow::qdaClassification, data = X)
                  parallel::stopCluster(cl)
                }
                qda.assignation <- qda.assignation[!is.na(qda.assignation)]
                qda.assignation <- matrix(unlist(qda.assignation), ncol = length(qda.assignation), byrow = FALSE)
                qda.groups <- Rfast::rowMaxs(qda.assignation)
                unique.qda.groups <- sort(unique(qda.groups))
                qda.groups <- factor(qda.groups)
                levels(qda.groups) <- unlist(lapply(assigned.training.cytometry, function(x) x$type))[!is.na(qda.assignation)]

                sys.time.step3.1 <- Sys.time()
                time.dif.3 <- difftime(sys.time.step3.1, sys.time.step3.0)
                print(paste("step 3:", time.dif.3, units(time.dif.3), sep = " "))
                t1 <- Sys.time()
                print(difftime(t1, t0))
                return(list(cluster = qda.groups, clusterings = tclust.results, assigned.template.index = c(assigned.template.index,
                  assigned.training.cytometry.index)))
            } else {
                assigned.training.cytometry <- database[cytos.cluster == assigned.template.index][[assigned.training.cytometry.index]]
                # ytest = as.factor(citometria$`Population ID (name)`)
                dim.cyto <- dim(assigned.training.cytometry)[2]
                y <- droplevels(as.factor(dplyr::pull(assigned.training.cytometry, dim.cyto)))
                training.rforest <- randomForest::randomForest(assigned.training.cytometry[, 1:(dim.cyto - 1)], y = y, keep.forest = TRUE)
                random.forest.groups <- predict(training.rforest, X)
                sys.time.step3.1 <- Sys.time()
                time.dif.3 <- difftime(sys.time.step3.1, sys.time.step3.0)
                print(paste("step 3:", time.dif.3, units(time.dif.3), sep = " "))
                t1 <- Sys.time()
                print(difftime(t1, t0))
                return(list(cluster = random.forest.groups, clusterings = tclust.results, assigned.template.index = c(assigned.template.index,
                  assigned.training.cytometry.index)))
            }
        }
    }
}
