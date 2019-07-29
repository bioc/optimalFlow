\name{estimClustPar}
\alias{estimClustPar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{estimClustPar
%%  ~~function to do ... ~~
}
\description{Obtain the best values for the model parameters, given data, input parameters and an assigment.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
estimClustPar(X, iter, pa)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{Points in a multidimiensional space.
%%     ~~Describe \code{X} here~~
}
  \item{iter}{Current solution obtained by tclust_.
%%     ~~Describe \code{iter} here~~
}
  \item{pa}{Parameters for using in tclust_.
%%     ~~Describe \code{pa} here~~
}
}
\value{An object used internally in tclust_ with the best model parameters.
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
## tclust_ is the function which obtain the clusters for tclust_H function

x=rbind(matrix(rnorm(100),ncol=2),matrix(rnorm(100)+2,ncol=2),matrix(rnorm(100)+4,ncol=2))
output3=tclust_( X=x , K=3 , alpha = 0.05 , niter = 20 , Ksteps=10 ,
equal.weights = FALSE, restr.cov.value = "eigen" ,
maxfact_e = 5 , zero.tol = 1e-16 ,  trace = 0 ,
sol_ini_p = FALSE ,   sol_ini=NA )

##restr.diffax <- function (iter, pa)
##	Apply constraints to covariance matrices
iter=output3$iter
pa=output3$pa

## estimClustPar  obtains the best values for the parameters,
##given data, input parameters and an assigment.
output4=estimClustPar (X=x, iter, pa)
output4$center
output4$sigma
output4$cw
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
