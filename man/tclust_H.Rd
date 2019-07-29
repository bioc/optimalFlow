\name{tclust_H}
\alias{tclust_H}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{tclust_H
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
tclust_H(x, k = 3, alpha = 0.05, nstart = 50, iter.max = 20, restr = "eigen", restr.fact = 12, sol_ini_p = FALSE, sol_ini = NA, equal.weights = FALSE, center = center, scale = scale, store.x = TRUE, drop.empty.clust = TRUE, trace = 0, warnings = 3, zero.tol = 1e-16)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{k}{
%%     ~~Describe \code{k} here~~
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
}
  \item{nstart}{
%%     ~~Describe \code{nstart} here~~
}
  \item{iter.max}{
%%     ~~Describe \code{iter.max} here~~
}
  \item{restr}{
%%     ~~Describe \code{restr} here~~
}
  \item{restr.fact}{
%%     ~~Describe \code{restr.fact} here~~
}
  \item{sol_ini_p}{
%%     ~~Describe \code{sol_ini_p} here~~
}
  \item{sol_ini}{
%%     ~~Describe \code{sol_ini} here~~
}
  \item{equal.weights}{
%%     ~~Describe \code{equal.weights} here~~
}
  \item{center}{
%%     ~~Describe \code{center} here~~
}
  \item{scale}{
%%     ~~Describe \code{scale} here~~
}
  \item{store.x}{
%%     ~~Describe \code{store.x} here~~
}
  \item{drop.empty.clust}{
%%     ~~Describe \code{drop.empty.clust} here~~
}
  \item{trace}{
%%     ~~Describe \code{trace} here~~
}
  \item{warnings}{
%%     ~~Describe \code{warnings} here~~
}
  \item{zero.tol}{
%%     ~~Describe \code{zero.tol} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, k = 3, alpha = 0.05, nstart = 50, iter.max = 20,
    restr = c("eigen", "deter", "sigma"), restr.fact = 12, sol_ini_p = FALSE,
    sol_ini = NA, equal.weights = FALSE, center = center, scale = scale,
    store.x = TRUE, drop.empty.clust = TRUE, trace = 0, warnings = 3,
    zero.tol = 1e-16)
{
    if (sol_ini_p == TRUE)
        if (is.numeric(sol_ini$cluster))
            sol_ini$assig = sol_ini$cluster
        else {
            sol_ini$cw = sol_ini$weights
            sol_ini$center = t(sol_ini$centers)
            sol_ini$sigma = sol_ini$cov
        }
    tc_ = tclust_(X = x, K = k, alpha = alpha, niter = nstart,
        Ksteps = iter.max, equal.weights = equal.weights, restr.cov.value = restr,
        maxfact_e = restr.fact, maxfact_d = Inf, m = 2, zero.tol = 1e-16,
        trace = 0, opt = "HARD", sol_ini_p = sol_ini_p, sol_ini = sol_ini,
        tk = FALSE)
    tc = list()
    tc$centers = t(tc_$iter$center)
    tc$cov = tc_$iter$sigma
    tc$cluster = tc_$iter$assig
    tc$par = tc_$pa
    tc$weights = tc_$iter$cw
    tc$obj = tc_$iter$obj
    return(tc)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
