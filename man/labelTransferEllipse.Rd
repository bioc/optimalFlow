\name{labelTransferEllipse}
\alias{labelTransferEllipse}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{labelTransferEllipse
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
labelTransferEllipse(i, test.cytometry.ellipses, training.cytometries.barycenter, equal.weights = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{i}{A dummy variable, should be any integral. Ment for use with lapply.}
  \item{test.cytometry.ellipses}{A test clustering viewed as a mixture of multivariate normal distributions.}
  \item{training.cytometries.barycenter}{ A training partition viewed as a mixture of multivariate normal distributions.}
  \item{equal.weights}{If True, weights assigned to every cluster in a partion are uniform (1/number of clusters) when calculating the similarity distance. If False, weights assigned to clusters are the proportions of points in every cluster compared to the total amount of points in the partition.}
%%     ~~Describe \code{x} here~~
}
\value{
A fuzzy relabeling consistent of a transportation plan.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\examples{
partition1 = list(list(mean = c(1,1), cov = diag(1,2), weight = 0.5, type = "1"), list(mean = c(-1,-1), cov = diag(1,2), weight = 0.5, type = "2"))
partition2 = list(list(mean = c(1,1), cov = diag(1,2), weight = 0.5, type = "a"), list(mean = c(-1,-1), cov = diag(1,2), weight = 0.5, type = "b"))
labelTransferEllipse(1, partition2, partition1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
