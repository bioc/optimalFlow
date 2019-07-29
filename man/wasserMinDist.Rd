\name{wasserMinDist}
\alias{wasserMinDist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{wasserMinDist
%%  ~~function to do ... ~~
}
\description{For two lists of multivariate normals calcualtes the closest member, in wasserstein distance, of the later list to each element of the former.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
wasserMinDist(points,centres)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{points}{List of multivariate normals, where each element is a list with values mean and cov.}
  \item{centres}{List of multivariate normals, where each element is a list with values mean and cov.}
%%     ~~Describe \code{x} here~~

}
\value{
A matrix where each column idicates the distance between the respective entry in points and the closest element in centers and the index of this closest element.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\examples{
normals = list(list(mean = c(1,1), cov = diag(2,2)),
list(mean = c(1,1), cov = diag(1,2)), list(mean = c(3,3), cov = diag(1,2)))
k_barycenter = kcenter(normals, 2, c(1,1,2))$kcenters
wasserMinDist(normals, k_barycenter)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
