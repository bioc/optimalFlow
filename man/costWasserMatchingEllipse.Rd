\name{costWasserMatchingEllipse}
\alias{costWasserMatchingEllipse}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{costWasserMatchingEllipse
%%  ~~function to do ... ~~
}
\description{Calculates similarity distance based on 2-Wassertein distance between mixtures of multivariate normal distributions.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
costWasserMatchingEllipse(test.cytometry, training.cytometries, equal.weights = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{test.cytometry}{A clusetering represented as a list of clusters. Each cluster is a list with elements mean, cov, weight and type.}
  \item{training.cytometries}{A list of clusterings with the same format as test.cytometry.}
  \item{equal.weights}{If True, weights assigned to every cluster in a partion are uniform (1/number of clusters) when calculating the similarity distance. If False, weights assigned to clusters are the proportions of points in every cluster compared to the total amount of points in the partition.}
%%     ~~Describe \code{x} here~~

}
\value{
 Returns a vector representing the similarity distance between test.cytometry and the elements in training.cytometries.
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
partition1 = list(list(mean = c(1,1), cov = diag(1,2), weight = 0.5, type = "1"),
list(mean = c(-1,-1), cov = diag(1,2), weight = 0.5, type = "2"))
partition2 = list(list(list(mean = c(1,-1), cov = diag(1,2),
weight = 0.5, type = "1"), list(mean = c(-1,1), cov = diag(1,2), weight = 0.5, type = "2")))
costWasserMatchingEllipse(partition1, partition2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
