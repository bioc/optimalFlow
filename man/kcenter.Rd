\name{kcenter}
\alias{kcenter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{kcenter
%%  ~~function to do ... ~~
}
\description{Calculates the k-barycenter of a list on multivaraite normals  for a given the number of clusters and an assignation of each normal to the respective group.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
kcenter(points, kk, center.asigned)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{points}{List of multivariate normals, where each element is a list with values mean and cov.
%%     ~~Describe \code{points} here~~
}
  \item{kk}{The number k of groups for the k-barycenter.
%%     ~~Describe \code{kk} here~~
}
  \item{center.asigned}{A vector indicating to which group each normal should be assigned.
%%     ~~Describe \code{center.asigned} here~~
}
}
\value{
\item{kcenters}{A list of the elements of the k-barycenter, i.e., a list of k lists containing means and covariances.}
\item{t.variation}{The trimmed wasserstein variation.}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\examples{
normals = list(list(mean = c(1,1), cov = diag(2,2)),
list(mean = c(1,1), cov = diag(1,2)), list(mean = c(3,3), cov = diag(1,2)))
kcenter(normals, 2, c(1,1,2))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
