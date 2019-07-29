\name{wassersteinKBarycenter}
\alias{wassersteinKBarycenter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{wassersteinKBarycenter
%%  ~~function to do ... ~~
}
\description{ A wrapper for calculating K-barycenters of multivariate normal distributions with the 2-Wasserstein distance.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
wassersteinKBarycenter(i = 1, k, alpha = 0, initialization = "rnd", pooled.clusters)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{i}{A dummy variable ment for use with apply.}
  \item{k}{Number k of elements in the k-barycenter.}
  \item{alpha}{Level of trimming.}
  \item{initialization}{Type of initialization in c("rnd", "plus-plus"). "rnd" makes the common random initilaization while "plus-plus" initializes in a similar fashion to k-means++.}
  \item{pooled.clusters}{List of multivariate normals for which the trimmed k-barycenter should be performed.}
%%     ~~Describe \code{x} here~~
}
\value{
\item{wasserstein.var}{Wasserstein variation.}
\item{wasserstein.k.barycenter}{List with three elements. Variacion_wasser is Waserstein variation. Baricentro is a list of k elements, each of which is a member of the k-barycenter. Cluster is the assignation to each barycenter of the original entries.}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
E del Barrio, H Inouzhe, JM Loubes, C Matran and A Mayo-Iscar. (2019) optimalFlow: Optimal-transport approach to flow cytometry gating and population matching. arXiv:1907.08006
}
\examples{
normals = list(list(mean = c(1,1), cov = diag(2,2)),
list(mean = c(1,1), cov = diag(1,2)), list(mean = c(3,3), cov = diag(1,2)))
wkb = wassersteinKBarycenter(1, 2, 0, "rnd", normals)
print(wkb$wasserstein.var)
print(wkb$wasserstein.k.barycente)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
