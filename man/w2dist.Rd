\name{w2dist}
\alias{w2dist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{w2dist
%%  ~~function to do ... ~~
}
\description{The 2-Wasserstein distance between two multivariate normal distributions
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
w2dist(P,Q)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{P}{A multivariate normal distribution given as a list with arguments mean and cov.}
  \item{Q}{A multivariate normal distribution given as a list with arguments mean and cov.}
%%     ~~Describe \code{x} here~~
}
\value{
The 2-Wasserstein distance between the two distributions.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\examples{
P = list(mean = c(1,1), cov = diag(1,2))
Q = list(mean = c(0,0), cov = 1.1*diag(1,2))
w2dist(P,Q)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
