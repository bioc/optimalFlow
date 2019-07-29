\name{tclustWithInitialization}
\alias{tclustWithInitialization}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{tclustWithInitialization
%%  ~~function to do ... ~~
}
\description{A wrapper for the function tclust_H.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
tclustWithInitialization(initialization, cytometry, i.sol.type = "points",
trimming = 0.05, restr.fact = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{initialization}{Initial solution for parameters provided by the user. Can be a matrix of data containing observations anc cluster assignations or can be a list spesifying a multivariate mixture of gaussians.}
  \item{cytometry}{A matrix or data.frame of dimension n x p, containing the observations (row-wise).}
  \item{i.sol.type}{Type of initial solutions in c("points", "barycenters"). "points" refers to a classified data matrix, while "barycenters" to a multivariate mixture.}
  \item{trimming}{The proportion of observations to be trimmed.}
  \item{restr.fact}{The constant restr.fact >= 1 constrains the allowed differences among group scatters. Larger values imply larger differences of group scatters, a value of 1 specifies the strongest restriction.}

%%     ~~Describe \code{x} here~~

}

\value{
\item{cluster}{A numerical vector of size n containing the cluster assignment for each observation. Cluster names are integer numbers from 1 to k, 0 indicates trimmed observations.}
\item{n_clus}{Number of clusters actually found.}
\item{obj}{he value of the objective function of the best (returned) solution.}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Fritz, H., Garcia-Escudero, L. A., & Mayo-Iscar, A. (2012). tclust: An r package for a trimming approach to cluster analysis. Journal of Statistical Software, 47(12), 1-26. E del Barrio, H Inouzhe, JM Loubes, C Matran and A Mayo-Iscar. (2019) optimalFlow: Optimal-transport approach to flow cytometry gating and population matching. arXiv:1907.08006.
}

\examples{
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

output2 = tclustWithInitialization(data.frame(x, output$cluster), x, "points", 0.05, 10)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
