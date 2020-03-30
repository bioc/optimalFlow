#' f1Score
#'
#' Calculates the F1 score fore each group in a partition.
#'
#' @param clustering The labels of the new classification.
#' @param cytometry Data of the clustering, where the last variable contains the original labels.
#' @param noise.cells An array of labels to be considered as noise.
#'
#' @return A matrix where the first row is the F1 score, the second row is the Precision and the third row is the Recall.
#'
#' @examples
#' f1Score(dplyr::pull(Cytometry3[c(sample(1:250,250),251:(dim(Cytometry3)[1])),],11),
#'         Cytometry3, noise.types)
#'
#' @import optimalFlowData
#' @import dplyr
#'
#' @references E del Barrio, H Inouzhe, JM Loubes, C Matran and A Mayo-Iscar. (2019) optimalFlow: Optimal-transport approach to flow cytometry gating and population matching. arXiv:1907.08006
#'
#' @export
#'
f1Score <- function(clustering, cytometry, noise.cells) {
    dim.cyto <- dim(cytometry)[2]
    prediction_qda_19_con_18 <- as.character(clustering)
    prediction_qda_19_con_18[prediction_qda_19_con_18 %in% noise.cells] <- rep("noise", sum(prediction_qda_19_con_18 %in% noise.cells))

    clustering <- prediction_qda_19_con_18
    lev_qda <- unique(clustering)
    p <- array(dim = c(1, length(lev_qda)))
    r <- p
    t <- 0
    for (i in lev_qda) {
        t <- t + 1
        if (i == "noise") {
            p[t] <- sum(clustering == i & dplyr::pull(cytometry, dim.cyto) %in% noise.cells)/sum(clustering == i)
            r[t] <- sum(clustering == i & dplyr::pull(cytometry, dim.cyto) %in% noise.cells)/sum(dplyr::pull(cytometry, dim.cyto) %in%
                noise.cells)
        } else {
            p[t] <- sum(clustering == i & dplyr::pull(cytometry, dim.cyto) %in% i)/sum(clustering == i)
            r[t] <- sum(clustering == i & dplyr::pull(cytometry, dim.cyto) %in% i)/sum(dplyr::pull(cytometry, dim.cyto) %in% i)
        }
    }
    F1_score <- 2 * p * r/(p + r)
    F1_score_tot <- rbind(F1_score, p, r)
    cels_original <- names(table(dplyr::pull(cytometry, dim.cyto)))[which(table(dplyr::pull(cytometry, dim.cyto)) > 0)]
    cel_no_encontradas <- cels_original[!(cels_original %in% prediction_qda_19_con_18) & !(cels_original %in% noise.cells)]
    F1_score_tot <- cbind(F1_score_tot, matrix(0, ncol = length(cel_no_encontradas), nrow = 3))
    colnames(F1_score_tot) <- c(lev_qda, cel_no_encontradas)
    rownames(F1_score_tot) <- c("F1-score", "Precision", "Recall")
    return(F1_score_tot)
}
