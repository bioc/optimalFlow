#' cytoPlotDatabase
#'
#' A plot wrapper for a database (list) of cytometries as a mixture of multivariate normals as used in optimalFlowTemplates.
#'
#' @param database.cytometries.as.mixtures A list where each component is a mixture distribution. That is, each component is a list, where each element contains the parameters of a component of the mixture as a list with entries: mean, cov, weight and type.
#' @param dimensions A vector containing the two variables on which to perform the projection.
#' @param xlim the x limits (x1, x2) of the plot. Note that x1 > x2 is allowed and leads to a ‘reversed axis’. The default value, NULL, indicates that the range of the finite values to be plotted should be used.
#' @param ylim the y limits of the plot.
#' @param xlab a label for the x axis, defaults to a description of x.
#' @param ylab a label for the y axis, defaults to a description of y.
#' @param colour If TRUE plots elements of a mixture distribution in different colours. If FALSE plots them in black.
#'
#' @return A two dimensional plot of ellipses containing the 95% of probability for the respectives components of the mixture distributions.
#' @examples
#' database <- buildDatabase(
#'   dataset_names = paste0('Cytometry', c(2:5, 7:9, 12:17, 19, 21)),
#'   population_ids = c('Monocytes', 'CD4+CD8-', 'Mature SIg Kappa', 'TCRgd-'))
#' templates.optimalFlow <-
#'   optimalFlowTemplates(
#'     database = database, templates.number = 5, cl.paral = 1
#'   )
#' cytoPlotDatabase(templates.optimalFlow$database.elliptical[which(templates.optimalFlow$clustering == 3)], dimensions = c(4,3), xlim = c(0, 8000), ylim = c(0, 8000), xlab = "", ylab = "")
#'
#' @import ellipse
#' @importFrom graphics lines
#' @export
cytoPlotDatabase <- function (database.cytometries.as.mixtures, dimensions = c(1,2), xlim = c(0,8000),
                      ylim = c(0,8000), xlab = "", ylab = "", colour = TRUE) {
  n.cytometries <- length(database.cytometries.as.mixtures)
  cytometry.as.mixture <- database.cytometries.as.mixtures[[1]]
  n.cells <- length(cytometry.as.mixture)
  plot(
    ellipse(
      cytometry.as.mixture[[1]]$cov[dimensions, dimensions],
      centre = cytometry.as.mixture[[1]]$mean[dimensions]
    ),
    xlim = xlim, ylim = ylim, col = 1, type = "l",
    xlab = xlab, ylab = ylab
  )
  if (colour) {
    for (j in 2:n.cells){
      lines(
        ellipse(
          cytometry.as.mixture[[j]]$cov[dimensions, dimensions],
          centre = cytometry.as.mixture[[j]]$mean[dimensions]
        ),
        col = j
      )
    }
    for (i in 2:n.cytometries){
      cytometry.as.mixture <- database.cytometries.as.mixtures[[i]]
      n.cells <- length(cytometry.as.mixture)
      for (j in 1:n.cells){
        lines(
          ellipse(
            cytometry.as.mixture[[j]]$cov[dimensions, dimensions],
            centre = cytometry.as.mixture[[j]]$mean[dimensions]
          ),
          col = j
        )
      }
    }
  } else {
    for (j in 2:n.cells){
      lines(
        ellipse(
          cytometry.as.mixture[[j]]$cov[dimensions, dimensions],
          centre = cytometry.as.mixture[[j]]$mean[dimensions]
        ),
        col = 1
      )
    }
    for (i in 2:n.cytometries){
      cytometry.as.mixture <- database.cytometries.as.mixtures[[i]]
      n.cells <- length(cytometry.as.mixture)
      for (j in 1:n.cells){
        lines(
          ellipse(
            cytometry.as.mixture[[j]]$cov[dimensions, dimensions],
            centre = cytometry.as.mixture[[j]]$mean[dimensions]
          ),
          col = 1
        )
      }
    }
  }
}
