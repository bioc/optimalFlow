#' voteLabelTransfer
#'
#' A wrapper for doing either labelTransfer or labelTransferEllipse.
#'
#' @param type 'points' indicates use of labelTransfer; 'ellipses' of labelTransferEllipse.
#' @param test.partition Only when type = 'points'. Labels of a partition of the test data.
#' @param test.cytometry Only when type = 'points'. Test data, a dataframe without labels.
#' @param test.partition.ellipse Only when type = 'ellipses'. A test clustering viewed as a mixture of multivariate normal distributions.
#' @param training.cytometries Only when type = 'points'. List of partitions, where each partition is a dataframe wher the last column contains the labels of the partition.
#' @param training.cytometries.barycenter Only when type = 'ellipses'. A training partition viewed as a mixture of multivariate normal distributions.
#' @param test Only when type = 'ellipses'. A dummy variable, should be any integral. Ment for use with lapply.
#' @param op.syst Type of system, takes values in c('unix', 'windows').
#' @param cl.paral Number of cores to be used in parallel procedures.
#' @param equal.weights If True, weights assigned to every cluster in a partion are uniform (1/number of clusters) when calculating the similarity distance. If False, weights assigned to clusters are the proportions of points in every cluster compared to the total amount of points in the partition.
#'
#' @return A list containing:
#' \describe{
#'  \item{final.vote}{A list for the votes on each cell.}
#'  \item{complete.vote}{A more complete list for the votes on each cell.}
#' }
#'
#' @examples
#' \donttest{
#' data.example <- data.frame(v1 = c(rnorm(50, 2, 1), rnorm(50, -2, 1)),
#'                           v2 = c(rnorm(50, 2, 1), rnorm(50, -2, 1)), id = c(rep(0, 50), rep(1, 50)))
#' test.labels <- c(rep('a', 50), rep('b', 50))
#' voteLabelTransfer(test.partition = test.labels, test.cytometry = data.example[, 1:2],
#'                   training.cytometries = list(data.example), op.syst = .Platform$OS.type)$final.vote[[1]]
#' }
#' @export
#'
voteLabelTransfer <- function(type = "points", test.partition, test.cytometry, test.partition.ellipse, training.cytometries, training.cytometries.barycenter,
    test = 1, op.syst, cl.paral = 1, equal.weights = FALSE) {
    if (type == "points") {
        transported.labels <- list()
        if (op.syst == "windows") {
            cl <- parallel::makePSOCKcluster(cl.paral)
            for (j in test) {
                t0 <- Sys.time()
                transported.labels[[as.character(j)]] <- parallel::parLapply(cl, training.cytometries, optimalFlow::labelTransfer,
                  test.cytometry, test.partition, equal.weights)
                t1 <- Sys.time()
                print(paste(t1 - t0, units(difftime(t1, t0))))
            }
            parallel::stopCluster(cl)
        } else {
            if (op.syst == "unix") {
                for (j in test) {
                  t0 <- Sys.time()
                  transported.labels[[as.character(j)]] <- parallel::mclapply(training.cytometries, optimalFlow::labelTransfer,
                    test.cytometry, test.partition, equal.weights, mc.cores = cl.paral)
                  t1 <- Sys.time()
                  print(paste(t1 - t0, units(difftime(t1, t0))))
                }

            }
        }
    } else {
        if (type == "ellipses") {
            test.partition <- 1:length(test.partition.ellipse)
            training.cytometries <- 1
            transported.labels <- list()
            for (j in test) {
                t0 <- Sys.time()
                transported.labels[[as.character(j)]] <- lapply(1, optimalFlow::labelTransferEllipse, test.partition.ellipse, training.cytometries.barycenter,
                  equal.weights)
                t1 <- Sys.time()
                print(paste(t1 - t0, units(difftime(t1, t0))))
            }
        }
    }

    vote.transport.labels <- list()
    complete.vote.transport.labels <- list()
    for (i in as.character(test)) {
        transp_lab_cito <- list()
        transp_lab_cito_complete <- list()
        for (test.cell in unique(transported.labels[[i]][[1]]$to)) {
            vote.on.test.cell <- c()
            for (j in 1:length(training.cytometries)) {
                receiving.cell.indexes <- which(transported.labels[[i]][[j]]$to == test.cell)
                sending.cells <- transported.labels[[i]][[j]]$from[receiving.cell.indexes]
                sending.cells.vote.proportion <- transported.labels[[i]][[j]]$mass[receiving.cell.indexes]/sum(transported.labels[[i]][[j]]$mass[receiving.cell.indexes])
                sending.cells.proportions <- array(dim = length(receiving.cell.indexes))
                t <- 0
                for (jj in receiving.cell.indexes) {
                  t <- t + 1
                  receiving.cell.indexes_jj <- which(transported.labels[[i]][[j]]$from == transported.labels[[i]][[j]]$from[jj])
                  sending.cells.proportions[t] <- transported.labels[[i]][[j]]$mass[jj]/sum(transported.labels[[i]][[j]]$mass[receiving.cell.indexes_jj])
                }
                vote.on.test.cell <- rbind(vote.on.test.cell, cbind(rep(j, t), sending.cells, sending.cells.vote.proportion, sending.cells.proportions))
            }
            vote.on.test.cell <- data.frame(id = as.integer(unlist(vote.on.test.cell[, 1])), cell = vote.on.test.cell[, 2], vote.proportion = as.double(unlist(vote.on.test.cell[,
                3])), original.proportion = as.double(unlist(vote.on.test.cell[, 4])))
            vote.on.test.cell$compound.proportion <- vote.on.test.cell$vote.proportion * vote.on.test.cell$original.proportion
            vote.on.test.cell.bis <- as.data.frame(vote.on.test.cell %>% group_by(cell) %>% summarise(compound.proportion = sum(compound.proportion),
                simple.proportion = sum(vote.proportion)))
            vote.on.test.cell.bis$compound.proportion <- vote.on.test.cell.bis$compound.proportion/sum(vote.on.test.cell.bis$compound.proportion)
            vote.on.test.cell.bis$simple.proportion <- vote.on.test.cell.bis$simple.proportion/sum(vote.on.test.cell.bis$simple.proportion)
            transp_lab_cito[[test.cell]] <- vote.on.test.cell.bis[order(-vote.on.test.cell.bis$compound.proportion), ]
            transp_lab_cito_complete[[test.cell]] <- vote.on.test.cell
        }
        vote.transport.labels[[i]] <- transp_lab_cito
        complete.vote.transport.labels[[i]] <- transp_lab_cito_complete
    }
    return(list(final.vote = vote.transport.labels, complete.vote = complete.vote.transport.labels))
}
