\name{GaussianBarycenters}
\alias{GaussianBarycenters}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{GaussianBarycenters
%%  ~~function to do ... ~~
}
\description{Finds the barricenter of a mixture of normals with the same mean.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GaussianBarycenters(matrices, weight)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{matrices}{A list of covariance matrices.
%%     ~~Describe \code{matrices} here~~
}
  \item{weight}{A vector of weights associated to the covariance matrices.
%%     ~~Describe \code{weight} here~~
}
}
\value{
\item{Barycenter}{Returns the barycenter, a covariance matrix.}
\item{Variation}{The Wasserstein Variation.}
\item{Num.iter}{THe number of iterations to achieve the stopping criteria.}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\examples{
GaussianBarycenters(list(diag(2,2),diag(1,2)), c(0.5,0.5))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
