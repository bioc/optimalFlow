\name{distGaussianCov}
\alias{distGaussianCov}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{distGaussianCov
%%  ~~function to do ... ~~
}
\description{Computes the componente relative to the covariances of the wasserstein distance.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
distGaussianCov(N1,N2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N1}{A multivariate normal distribution as list with elements mean and cov.}
  \item{N2}{A multivariate normal distribution as list with elements mean and cov.}
%%     ~~Describe \code{x} here~~

}
\value{
Wasserstein squared distance between normals with the same mean.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\examples{
distGaussianCov(list(mean = c(1,1), cov = diag(2,2)), list(mean = c(1,1), cov = diag(1,2)))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
