\name{voteLabelTransfer}
\alias{voteLabelTransfer}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{voteLabelTransfer
%%  ~~function to do ... ~~
}
\description{A wrapper for doing either labelTransfer or labelTransferEllipse
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{voteLabelTransfer(type = "points", test.partition, test.cytometry,
test.partition.ellipse, training.cytometries,training.cytometries.barycenter,
test = 1, op.syst, cl.paral = 1, equal.weights = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{type}{"points" indicates use of labelTransfer; "ellipses" of labelTransferEllipse.}
  \item{test.partition}{Only when type = "points". Labels of a partition of the test data.}
  \item{test.cytometry}{Only when type = "points". Test data, a dataframe without labels.}
  \item{test.partition.ellipse}{Only when type = "ellipses". A test clustering viewed as a mixture of multivariate normal distributions.}
  \item{training.cytometries}{Only when type = "points". List of partitions, where each partition is a dataframe wher the last column contains the labels of the partition.}
  \item{training.cytometries.barycenter}{Only when type = "ellipses". A training partition viewed as a mixture of multivariate normal distributions.}
  \item{test}{Only when type = "ellipses". A dummy variable, should be any integral. Ment for use with lapply.}
  \item{op.syst}{Type of system, takes values in c("unix", "windows").}
  \item{cl.paral}{Number of cores to be used in parallel procedures.}
  \item{equal.weights}{If True, weights assigned to every cluster in a partion are uniform (1/number of clusters) when calculating the similarity distance. If False, weights assigned to clusters are the proportions of points in every cluster compared to the total amount of points in the partition.}
%%     ~~Describe \code{x} here~~
}
\value{
\item{final.vote}{A list for the votes on each cell.}
\item{complete.vote}{A more complete list for the votes on each cell.}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\examples{
data.example = data.frame(v1 = c(rnorm(50,2,1), rnorm(50,-2,1)),
v2 = c(rnorm(50,2,1), rnorm(50,-2,1)), id = c(rep(0,50), rep(1,50)))
test.labels = c(rep("a",50), rep("b", 50))
voteLabelTransfer(test.partition = test.labels, test.cytometry = data.example[,1:2],
training.cytometries = list(data.example), op.syst = .Platform$OS.type)$final.vote[[1]]
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
