\name{estimationCellBarycenter}
\alias{estimationCellBarycenter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{estimationCellBarycenter
%%  ~~function to do ... ~~
}
\description{Estimates a Wasserstein barycenter for a cluster type using a collection of partitions.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
estimationCellBarycenter(cell, cytometries)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cell}{Name of the cluster of interes}
  \item{cytometries}{List of clusterings}
%%     ~~Describe \code{x} here~~

}
\value{
\item{mean}{Mean of the barycenter.}
\item{cov}{Covariance of the barycenter.}
\item{weight}{Weight associated to the barycenter.}
\item{type}{Type of the cluster.}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\examples{
partition1 = list(list(mean = c(1,1), cov = diag(1,2), weight = 0.5, type = "1"),
list(mean = c(-1,-1), cov = diag(1,2), weight = 0.5, type = "2"))
partition2 = list(list(mean = c(1,-1), cov = diag(1,2), weight = 0.5, type = "1"),
list(mean = c(-1,1), cov = diag(1,2), weight = 0.5, type = "2"))
cytometries = list(partition1, partition2)
estimationCellBarycenter("1",cytometries)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
