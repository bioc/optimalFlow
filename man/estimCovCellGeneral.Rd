\name{estimCovCellGeneral}
\alias{estimCovCellGeneral}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{estimCovCellGeneral
%%  ~~function to do ... ~~
}
\description{Estimation of mean and covariance for a label in a partition.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
estimCovCellGeneral(cell, cytometry, labels, type = "standard", alpha = 0.85)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cell}{Labell of the clsuter of interest.}
  \item{cytometry}{Data of the partition, wthout labels.}
  \item{labels}{Labels of the partition.}
  \item{type}{How to estimate covariance matrices of a cluster. "standard" is for using cov(), while "robust" is for using robustbase::covMcd.}
  \item{alpha}{Only when type = "robust". Indicates the value of alpha in robustbase::covMcd.}
%%     ~~Describe \code{x} here~~

}
\value{
\item{mean}{Mean of the cluster.}
\item{cov}{Covariance of the cluster.}
\item{weight}{Weight associated to the cluster.}
\item{type}{Type of the cluster.}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\examples{
estimCovCellGeneral("Basophils", Cytometry1[,1:10], Cytometry1[,11])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
