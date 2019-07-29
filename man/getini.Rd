\name{getini}
\alias{getini}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{getini
%%  ~~function to do ... ~~
}
\description{Calculates the initial cluster sizes for initializing tclust_.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getini(K, no.trim)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{K}{ Number of groups.
%%     ~~Describe \code{K} here~~
}
  \item{no.trim}{Indicates if trimming is allowed.
%%     ~~Describe \code{no.trim} here~~
}
}
\value{
A vector of length K indicating indices of the points to be considered for initialization.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}

\examples{
## gives a random vector from a K dimensional multinomial(no.trim, pi.ini)
## with pi.ini a random vector build with random values from uniform (0,1)
v=getini  (K=3, no.trim=100)
v
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
