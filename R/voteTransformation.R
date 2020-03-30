#' voteTransformation
#'
#' Transforming votes obtained by using optimalFlowTemplates + OptimalFlowClassification with consenus.method in c('hierarchical', 'k-barycenter') and classif.method = 'matching' and cost.function = 'ellipses' to an appropriate format for using f1ScoreVoting.
#'
#' @param vote.0 Values obtained by voteLabelTransfer.
#' @param vote.1 Original proportions of the clusters after the template obtention.
#'
#' @return A list for the votes on each cell.
#'
#' @examples
#' vote.0 <- list('1' = data.frame(cell = c(1, 2), 'compound.proportion' = c(0.7, 0.3),
#'                                'simple.proportion'= c(0.7, 0.3)), '2' = data.frame(cell = c(1, 2),
#'                                                                                   'compound.proportion' = c(0.3, 0.7), 'simple.proportion'= c(0.3, 0.7)))
#' vote.1.1 <- t(c(0.8, 0.2))
#' names(vote.1.1) <- c('A', 'B')
#' vote.1.2 <- t(c(0.2, 0.8))
#' names(vote.1.2) <- c('A', 'B')
#' vote.1 <- list(vote.1.1, vote.1.2)
#' voteTransformation(vote.0, vote.1)
#'
#' @noRd
#'
voteTransformation <- function(vote.0, vote.1) {
    final.vote <- lapply(names(vote.0), function(i) {
        cells.0 <- as.character(vote.0[[i]]$cell)
        cells.1 <- vote.1[as.integer(cells.0)]
        joint.data.frame <- data.frame()
        joint.data.frame.1 <- data.frame()
        for (j in 1:length(cells.1)) {
            temp.data.frame <- data.frame(cell = names(cells.1[[j]]), simple.proportion = vote.0[[i]]$simple.proportion[j] * as.vector(cells.1[[j]]), 
                compound.proportion = vote.0[[i]]$compound.proportion[j] * as.vector(cells.1[[j]]))
            joint.data.frame = rbind(joint.data.frame, temp.data.frame)
        }
        cell.vote <- joint.data.frame %>% group_by(cell) %>% summarise(simple.proportion = sum(simple.proportion), compound.proportion = sum(compound.proportion))
        cell.vote$simple.proportion <- cell.vote$simple.proportion/sum(cell.vote$simple.proportion)
        cell.vote$compound.proportion <- cell.vote$compound.proportion/sum(cell.vote$compound.proportion)
        cell.vote <- as.data.frame(dplyr::arrange(cell.vote, desc(compound.proportion)))
    })
    names(final.vote) <- names(vote.0)
    return(final.vote)
}
