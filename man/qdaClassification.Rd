\name{qdaClassification}
\alias{qdaClassification}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{qdaClassification
%%  ~~function to do ... ~~
}
\description{ GIves quadratic discriminant scores to the poiints in data for a multivariate normal.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
qdaClassification(normal, data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{normal}{A list with arguments mean, covaruance and weight.}
  \item{data}{Data frame or matrix on which to perform qda.}
%%     ~~Describe \code{x} here~~

}
\value{
A score for each point.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\examples{
data.qda = cbind(rnorm(50), rnorm(50))
exp(qdaClassification(list(mean = c(0,0), cov = diag(1,2), weight = 1), data.qda))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
