\name{ssclmat}
\alias{ssclmat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ssclmat
%%  ~~function to do ... ~~
}
\description{Extract a matrix from the object containing covariance matrices. Ment for use inside tclust_ function.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ssclmat(x, k)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An object used inside tclust_ containing covariance matrices.
%%     ~~Describe \code{x} here~~
}
  \item{k}{An integer value giving a location in an array.
%%     ~~Describe \code{k} here~~
}
}
\value{
A covariance matrix.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\examples{
######EXAMPLE tclust_
##### tclust_ is the function which obtain the clusters for tclust_H function

x=rbind(matrix(rnorm(100),ncol=2),matrix(rnorm(100)+2,ncol=2),matrix(rnorm(100)+4,ncol=2))
output3=tclust_( X=x , K=3 , alpha = 0.05 , niter = 20 , Ksteps=10 ,
equal.weights = FALSE, restr.cov.value = "eigen" ,
maxfact_e = 5 , zero.tol = 1e-16 ,  trace = 0 ,
sol_ini_p = FALSE ,   sol_ini=NA )

##restr.diffax <- function (iter, pa)
##	Apply constraints to covariance matrices
iter=output3$iter
pa=output3$pa
pa$maxfact_e=1.1
iter2=restr.diffax (iter, pa)

##EXAMPLE	extract matrix from the object containing covariance matrices
##sclmat <- function (x, k) as.matrix (x[,,k])
ssclmat(iter2$sigma,k=1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
