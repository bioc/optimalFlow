\name{dmnorm}
\alias{dmnorm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{dmnorm
%%  ~~function to do ... ~~
}
\description{ The multivariate normal density.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
dmnorm(X, mu, sigma)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{ Multivariate points.
%%     ~~Describe \code{X} here~~
}
  \item{mu}{Mean.
%%     ~~Describe \code{mu} here~~
}
  \item{sigma}{COvariance matrix.
%%     ~~Describe \code{sigma} here~~
}
}
\value{
The value of the multivariate normal with mean mu and covariance sigma at the points in X.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}

\examples{
##	Multivariate normal density
##	Gives Multivariate normal density values for given  mu  and sigma
x=rbind(matrix(rnorm(100),ncol=2),matrix(rnorm(100)+2,ncol=2),matrix(rnorm(100)+4,ncol=2))
dmnorm (X=x,mu=c(0,0),sigma=diag(2))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
