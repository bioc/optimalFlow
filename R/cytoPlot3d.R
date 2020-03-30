#' cytoPlot3d
#'
#' A rgl::plot3d wrapper for cytometries as a mixture of multivariate normals as used in optimalFlowTemplates.
#'
#' @param cytometry.as.mixture A list, where each element contains the parameters of a component of the mixture as a list with entries: mean, cov, weight and type.
#' @param dimensions A vector containing the three variables on which to perform the projection.
#' @param xlim the x limits (x1, x2) of the plot. Note that x1 > x2 is allowed and leads to a ‘reversed axis’. The default value, NULL, indicates that the range of the finite values to be plotted should be used.
#' @param ylim the y limits of the plot.
#' @param zlim the z limits of the plot.
#' @param xlab a label for the x axis, defaults to a description of x.
#' @param ylab a label for the y axis, defaults to a description of y.
#' @param zlab a label for the z axis, defaults to a description of z.
#'
#' @return A three dimensional plot of ellipsoids containing the 95% of probability for the respectives components of the mixture distribution.
#' @examples
#' database <- buildDatabase(
#'   dataset_names = paste0('Cytometry', c(2:5, 7:9, 12:17, 19, 21)),
#'   population_ids = c('Monocytes', 'CD4+CD8-', 'Mature SIg Kappa', 'TCRgd-'))
#' templates.optimalFlow <-
#'   optimalFlowTemplates(
#'     database = database, templates.number = 5, cl.paral = 1
#'   )
#' # # To execute requires an actual monitor since it uses rgl.
#' # cytoPlot3d(templates.optimalFlow$templates[[3]], dimensions = c(4, 3, 9), xlim = c(0, 8000), ylim = c(0, 8000), zlim = c(0, 8000), xlab = "", ylab = "", zlab = ")
#'
#' @import rgl
#' @export
cytoPlot3d <- function (cytometry.as.mixture, dimensions = c(1, 2), xlim = NULL,
                        ylim = NULL, zlim = NULL, xlab = NULL, ylab = NULL, zlab = NULL) {
  n.cells <- length(cytometry.as.mixture)
  rgl::plot3d(rgl::ellipse3d(cytometry.as.mixture[[1]]$
                               cov[dimensions, dimensions],
                             centre = cytometry.as.mixture[[1]]$
                               mean[dimensions]),
              xlim = xlim, ylim = ylim, zlim = zlim, alpha = 0.5,
              col = 1, xlab = xlab, ylab = ylab,
              zlab = zlab)
  for (j in 2:n.cells){
    rgl::plot3d(rgl::ellipse3d(cytometry.as.mixture[[j]]$cov[dimensions, dimensions],
                               centre = cytometry.as.mixture[[j]]$
                                 mean[dimensions]), alpha = 0.5, add = TRUE, col = j)
  }
}
