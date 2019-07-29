\name{tclust_}
\alias{tclust_}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{tclust_
%%  ~~function to do ... ~~
}
\description{Function that performs robust non spherical clustering, tclust, where initial solutions are allowed.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
tclust_(X, K, alpha = 0.05, niter = 20, Ksteps = 10, equal.weights = FALSE,
restr.cov.value = "eigen", maxfact_e = 5, zero.tol = 1e-16, trace = 0,
sol_ini_p = FALSE, sol_ini = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{Points in a multidimiensional space.
%%     ~~Describe \code{X} here~~
}
  \item{K}{Number of clusters.
%%     ~~Describe \code{K} here~~
}
  \item{alpha}{Level of trimming.
%%     ~~Describe \code{alpha} here~~
}
  \item{niter}{ Maximum number of iterations.
%%     ~~Describe \code{niter} here~~
}
  \item{Ksteps}{ Maximum number of K-steps.
%%     ~~Describe \code{Ksteps} here~~
}
  \item{equal.weights}{If all clusters should have equal weight in the objective function.
%%     ~~Describe \code{equal.weights} here~~
}
  \item{restr.cov.value}{Type of restriction on covariance matrices. Only "eigen" will be allowed.
%%     ~~Describe \code{restr.cov.value} here~~
}
  \item{maxfact_e}{Level determinant constraints.
%%     ~~Describe \code{maxfact_e} here~~
}
  \item{zero.tol}{The zero tolerance used. By default set to 1e-16..
%%     ~~Describe \code{zero.tol} here~~
}
  \item{trace}{Defines the tracing level, which is set to 0 by default. Tracing level 2 gives additional information on the iteratively decreasing objective function's value.
%%     ~~Describe \code{trace} here~~
}
  \item{sol_ini_p}{Initial solution for parameters provided by the user TRUE/FALSE, if TRUE is stored in sol_ini.
%%     ~~Describe \code{sol_ini_p} here~~
}
  \item{sol_ini}{Initial solution for parameters provided by the user.
%%     ~~Describe \code{sol_ini} here~~
}
}
\value{
\item{iter}{Contains the solutions of the clustering provided by tclust_.}
\item{pa}{Parametaers used in tslucs_.}
\item{X}{The input variables.}
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
x=rbind(matrix(rnorm(100),ncol=2),matrix(rnorm(100)+2,ncol=2),matrix(rnorm(100)+4,ncol=2))
output3=tclust_( X=x , K=3 , alpha = 0.05 , niter = 20 , Ksteps=10 ,
equal.weights = FALSE, restr.cov.value = "eigen" ,
maxfact_e = 5 , zero.tol = 1e-16 ,  trace = 0 ,
sol_ini_p = FALSE ,   sol_ini=NA )
## cluster assigment
output3$iter$assig
plot(x,col=output3$iter$assig+1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
