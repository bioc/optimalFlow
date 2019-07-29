\name{trimmedKBarycenter}
\alias{trimmedKBarycenter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{trimmedKBarycenter
%%  ~~function to do ... ~~
}
\description{Calculates K-barycenters of multivariate normal distributions with the 2-Wasserstein distance
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
trimmedKBarycenter(k, alpha0, type.ini = "rnd", reps.list)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{k}{Number k of elements in the k-barycenter.}
  \item{alpha0}{Level of trimming.}
  \item{type.ini}{Type of initialization in c("rnd", "plus-plus"). "rnd" makes the common random initilaization while "plus-plus" initializes in a similar fashion to k-means++.}
  \item{reps.list}{List of multivariate normals for which the trimmed k-barycenter should be performed.}
%%     ~~Describe \code{x} here~~

}
\value{
\item{variacion_wasser}{Waserstein variation.}
\item{baricentro}{A list of k elements, each of which is a member of the k-barycenter}
\item{cluster}{The assignment of the original entries to each member of the k-barycenter.}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\examples{
normals = list(list(mean = c(1,1), cov = diag(2,2)), list(mean = c(1,1),
cov = diag(1,2)), list(mean = c(3,3), cov = diag(1,2)))
trimmedKBarycenter(2, 0, "rnd", normals)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
