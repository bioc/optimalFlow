\name{tclust_H}
\alias{tclust_H}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{tclust_H
%%  ~~function to do ... ~~
}
\description{ A wrapper for the fucntion tclust_. Performs robust non spherical clustering, tclust, where initial solutions are allowed.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
tclust_H(x, k = 3, alpha = 0.05, nstart = 50, iter.max = 20, restr = "eigen",
restr.fact = 12, sol_ini_p = FALSE, sol_ini = NA, equal.weights = FALSE,
center = center, scale = scale, store.x = TRUE, drop.empty.clust = TRUE,
trace = 0, warnings = 3, zero.tol = 1e-16)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A matrix or data.frame of dimension n x p, containing the observations (row-wise).
%%     ~~Describe \code{x} here~~
}
  \item{k}{The number of clusters initially searched for.
%%     ~~Describe \code{k} here~~
}
  \item{alpha}{The proportion of observations to be trimmed.
%%     ~~Describe \code{alpha} here~~
}
  \item{nstart}{The number of random initializations to be performed. Only when sol_ini_p = FALSE.
%%     ~~Describe \code{nstart} here~~
}
  \item{iter.max}{The maximum number of concentration steps to be performed. The concentration steps are stopped, whenever two consecutive steps lead to the same data partition.
%%     ~~Describe \code{iter.max} here~~
}
  \item{restr}{The type of restriction to be applied on the cluster scatter matrices. Valid values are "eigen" (default).
%%     ~~Describe \code{restr} here~~
}
  \item{restr.fact}{The constant restr.fact >= 1 constrains the allowed differences among group scatters. Larger values imply larger differences of group scatters, a value of 1 specifies the strongest restriction.
%%     ~~Describe \code{restr.fact} here~~
}
  \item{sol_ini_p}{Initial solution for parameters provided by the user TRUE/FALSE, if TRUE is stored in sol_ini.
%%     ~~Describe \code{sol_ini_p} here~~
}
  \item{sol_ini}{Initial solution for parameters provided by the user.
%%     ~~Describe \code{sol_ini} here~~
}
  \item{equal.weights}{A logical value, specifying whether equal cluster weights (TRUE) or not (FALSE) shall be considered in the concentration and assignment steps.
%%     ~~Describe \code{equal.weights} here~~
}
  \item{center}{A center vector of length p which can optionally be specified for centering x before calculation.
%%     ~~Describe \code{center} here~~
}
  \item{scale}{A scale vector of length p which can optionally be specified for scaling x before calculation.
%%     ~~Describe \code{scale} here~~
}
  \item{store.x}{A logical value, specifying whether the data matrix x shall be included in the result structure. By default this value is set to TRUE.
%%     ~~Describe \code{store.x} here~~
}
  \item{drop.empty.clust}{Logical value specifying, whether empty clusters shall be omitted in the resulting object. (The result structure does not contain center and covariance estimates of empty clusters anymore. Cluster names are reassigned such that the first l clusters (l <= k) always have at least one observation.
%%     ~~Describe \code{drop.empty.clust} here~~
}
  \item{trace}{Defines the tracing level, which is set to 0 by default. Tracing level 2 gives additional information on the iteratively decreasing objective function's value.
%%     ~~Describe \code{trace} here~~
}
  \item{warnings}{The warning level (0: no warnings; 1: warnings on unexpected behavior; 2: warnings if restr.fact causes artificially restricted results).
%%     ~~Describe \code{warnings} here~~
}
  \item{zero.tol}{The zero tolerance used. By default set to 1e-16.
%%     ~~Describe \code{zero.tol} here~~
}
}
\details{
This iterative algorithm initializes k clusters randomly and performs "concentration steps" in order to improve the current cluster assignment. The number of maximum concentration steps to be performed is given by iter.max. For approximately obtaining the global optimum, the system is initialized nstart times and concentration steps are performed until convergence or iter.max is reached. When processing more complex data sets higher values of nstart and iter.max have to be specified (obviously implying extra computation time). However, if more then half of the iterations would not converge, a warning message is issued, indicating that nstart has to be increased.

The parameter restr defines the cluster's shape restrictions, which are applied on all clusters during each iteration. Options "eigen"/"deter" restrict the ratio between the maximum and minimum eigenvalue/determinant of all cluster's covariance structures to parameter restr.fact. Setting restr.fact to 1, yields the strongest restriction, forcing all eigenvalues/determinants to be equal and so the method looks for similarly scattered (respectively spherical) clusters. Option "sigma" is a simpler restriction, which averages the covariance structures during each iteration (weighted by cluster sizes) in order to get similar (equal) cluster scatters.
}
\value{
\item{centers}{A matrix of size p x k containing the centers (column-wise) of each cluster.}
\item{cov}{An array of size p x p x k containing the covariance matrices of each cluster.}
\item{cluster}{A numerical vector of size n containing the cluster assignment for each observation. Cluster names are integer numbers from 1 to k, 0 indicates trimmed observations.}
\item{par}{A list, containing the parameters the algorithm has been called with (x, if not suppressed by store.x = FALSE, k, alpha, restr.fact, nstart, KStep, and equal.weights).}
\item{weights}{A numerical vector of length k, containing the weights of each cluster.}
\item{obj}{he value of the objective function of the best (returned) solution.}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Fritz, H., Garcia-Escudero, L. A., & Mayo-Iscar, A. (2012). tclust: An r package for a trimming approach to cluster analysis. Journal of Statistical Software, 47(12), 1-26.
}
\examples{

## tclust_H if the function which gives clusters to the user.
## The main role of this function is to be an interface with
## the user using labels for the parameters similar to tclust
## function in tclust package

x=rbind(matrix(rnorm(100),ncol=2),matrix(rnorm(100)+2,ncol=2),
matrix(rnorm(100)+4,ncol=2))

## robust cluster obtention from a sample x asking for 3 clusters,
## trimming level 0.05 and constrain level 12

k=3;alpha = 0.05;restr.fact = 12
output=tclust_H (x=x, k = k, alpha = alpha, nstart = 50, iter.max = 20,
restr = "eigen", restr.fact = restr.fact, sol_ini_p = FALSE, sol_ini = NA,
equal.weights = FALSE, center = center, scale = scale, store.x = TRUE,
drop.empty.clust = TRUE, trace = 0, warnings = 3, zero.tol = 1e-16)

## cluster assigment

output$cluster
plot(x,col=output$cluster)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
