\name{labelTransfer}
\alias{labelTransfer}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{labelTransfer
%%  ~~function to do ... ~~
}
\description{Label transfer between a test partition and a training set of partitions.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
labelTransfer(training.cytometry, test.cytometry, test.partition, equal.weights = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{training.cytometry}{List of partitions, where each partition is a dataframe wher the last column contains the labels of the partition.}
  \item{test.cytometry}{Test data, a dataframe without labels.}
  \item{test.partition}{Labels of a partition of the test data.}
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
E del Barrio, H Inouzhe, JM Loubes, C Matran and A Mayo-Iscar. (2019) optimalFlow: Optimal-transport approach to flow cytometry gating and population matching. arXiv:1907.08006
}
\examples{
data.example = data.frame(v1 = c(rnorm(50,2,1), rnorm(50,-2,1)),
v2 = c(rnorm(50,2,1), rnorm(50,-2,1)), id = c(rep(0,50), rep(1,50)))
test.labels = c(rep("a",50), rep("b", 50))
labelTransfer(data.example, data.example[,1:2], test.labels)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
