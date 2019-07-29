\name{restr2_eigenv}
\alias{restr2_eigenv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{restr2_eigenv
%%  ~~function to do ... ~~
}
\description{Function for applying eigenvalue constraints. These are the typical constraints used in tclust.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
restr2_eigenv(autovalues, ni.ini, factor_e, zero.tol)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{autovalues}{Matrix containin eigenvalues.
%%     ~~Describe \code{autovalues} here~~
}
  \item{ni.ini}{Current sample size of the clusters.
%%     ~~Describe \code{ni.ini} here~~
}
  \item{factor_e}{The level of the constraints.
%%     ~~Describe \code{factor_e} here~~
}
  \item{zero.tol}{Toletance level.
%%     ~~Describe \code{zero.tol} here~~
}
}
\value{
The restricted eigenvalues
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Fritz, H., Garcia-Escudero, L. A., & Mayo-Iscar, A. (2012). tclust: An r package for a trimming approach to cluster analysis. Journal of Statistical Software, 47(12), 1-26.
}
\examples{
#restr2_eigenv <- function(autovalues, ni.ini, factor_e, zero.tol)
#gives optimal constrained eigenvalues
autovalues=matrix(c(2,3,4,1,2,3),nrow=2)
ni.ini=c(2,2,3)
factor_e=1.1
zero.tol=1e-9
autovalues_const= restr2_eigenv (autovalues, ni.ini, factor_e, zero.tol)
autovalues_const
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
