\name{TreatSingularity}
\alias{TreatSingularity}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{TreatSingularity
%%  ~~function to do ... ~~
}
\description{Shows a warning message when a singularity in the curent solution of tclust_ is found.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
TreatSingularity(iter, pa)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{iter}{Current solution obtained by tclust_
%%     ~~Describe \code{iter} here~~
}
  \item{pa}{Parameters for using in tclust_.
%%     ~~Describe \code{pa} here~~
}
}
\value{
A warning message.
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

##### It shows a warning message
#####	warning ("points in the data set are concentrated in k points after trimming ")
TreatSingularity (iter, pa)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
