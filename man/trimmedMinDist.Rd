\name{trimmedMinDist}
\alias{trimmedMinDist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{trimmedMinDist
%%  ~~function to do ... ~~
}
\description{For two lists of multivaraite normals, points and centers, returns the index of which element in centers is closest to each element in points when some trimming is allowed.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
trimmedMinDist(points, centres, alpha = 0.1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{points}{List of multivariate normals, where each element is a list with values mean and cov.}
  \item{centres}{List of multivariate normals, where each element is a list with values mean and cov.}
  \item{alpha}{Level of triming}
%%     ~~Describe \code{x} here~~

}
\value{
A vector with the index of which element in centers is closest to each element in points.
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
trimmedMinDist(normals,k_barycenter, 0)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
