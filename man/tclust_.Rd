\name{tclust_}
\alias{tclust_}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{tclust_
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
tclust_(X, K, alpha = 0.05, niter = 20, Ksteps = 10, equal.weights = FALSE, restr.cov.value = "eigen", maxfact_e = 5, zero.tol = 1e-16, trace = 0, sol_ini_p = FALSE, sol_ini = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{K}{
%%     ~~Describe \code{K} here~~
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
}
  \item{niter}{
%%     ~~Describe \code{niter} here~~
}
  \item{Ksteps}{
%%     ~~Describe \code{Ksteps} here~~
}
  \item{equal.weights}{
%%     ~~Describe \code{equal.weights} here~~
}
  \item{restr.cov.value}{
%%     ~~Describe \code{restr.cov.value} here~~
}
  \item{maxfact_e}{
%%     ~~Describe \code{maxfact_e} here~~
}
  \item{zero.tol}{
%%     ~~Describe \code{zero.tol} here~~
}
  \item{trace}{
%%     ~~Describe \code{trace} here~~
}
  \item{sol_ini_p}{
%%     ~~Describe \code{sol_ini_p} here~~
}
  \item{sol_ini}{
%%     ~~Describe \code{sol_ini} here~~
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
function (X, K, alpha = 0.05, niter = 20, Ksteps = 10, equal.weights = FALSE,
    restr.cov.value = "eigen", maxfact_e = 5, maxfact_d = 5,
    m = 2, zero.tol = 1e-16, trace = 0, opt = "HARD", sol_ini_p = FALSE,
    sol_ini = NA, tk = FALSE)
{
    if (!is.numeric(X))
        stop("parameter x: numeric matrix/vector expected")
    if (!is.matrix(X))
        X <- matrix(X, ncol = 1)
    n <- nrow(X)
    p <- ncol(X)
    no.trim <- floor(n * (1 - alpha))
    if (restr.cov.value == "sigma") {
        f.restr = restr.avgcov
        restr.deter = FALSE
    }
    if (p == 1) {
        f.restr = restr.diffax
        restr.deter = FALSE
    }
    if ((p > 1) & (restr.cov.value == "eigen")) {
        f.restr = restr.diffax
        restr.deter = FALSE
    }
    if ((p > 1) & (restr.cov.value == "deter")) {
        f.restr = restr.diffax
        restr.deter = TRUE
    }
    pa <- list(n = n, p = p, alpha = alpha, trimm = n - no.trim,
        no.trim = no.trim, K = K, equal.weights = equal.weights,
        zero.tol = zero.tol, trace = trace, m = m, restr.deter = restr.deter,
        maxfact_d = maxfact_d, maxfact_e = maxfact_e, opt = opt,
        sol_ini_p = sol_ini_p, sol_ini = sol_ini, Ksteps = Ksteps,
        niter = niter, tk = tk)
    iter <- list(obj = NA, assig = array(0, n), csize = array(NA,
        K), cw = rep(NA, K), sigma = array(NA, c(p, p, K)), center = array(NA,
        c(K, p)), code = NA, z_ij = matrix(0, nrow = n, ncol = K),
        lambda = rep(NA, K))
    best.iter <- list(obj = -Inf)
    if (pa$sol_ini_p == TRUE)
        pa$niter = 1
    for (j in 1:pa$niter) {
        if (pa$sol_ini_p == TRUE) {
            if (is.numeric(sol_ini$assig)) {
                iter$assig = sol_ini$assig
                iter$z_ij = 0 * iter$z_ij
                for (h in 1:pa$K) iter$z_ij[iter$assig == h,
                  h] = 1
                iter$csize <- tabulate(iter$assig, pa$K)
                iter$cw <- iter$csize/sum(iter$csize)
                iter <- estimClustPar(X, iter, pa)
                iter$center[is.na(iter$center)] = 0
            }
            else {
                iter$cw = sol_ini$cw
                iter$center = sol_ini$center
                iter$sigma = sol_ini$sigma
                iter$csize <- iter$cw * pa$n * (1 - pa$alpha)
            }
        }
        else iter <- InitClusters(X, iter, pa)
        iter$center[is.na(iter$center)] = 0
        for (i in 0:pa$Ksteps) {
            iter <- f.restr(iter = iter, pa = pa)
            if (iter$code == 0) {
                if (i > 0)
                  return(TreatSingularity(calcobj(X, iter, pa),
                    pa))
                else iter$sigma[, , ] = diag(pa$p)
            }
            if (pa$opt == "FUZZY")
                iter <- findClustAssig_f(X, iter, pa)
            else iter <- findClustAssig(X, iter, pa)
            if ((iter$code == 2) || (i == pa$Ksteps))
                break
            iter <- estimClustPar(X, iter, pa)
        }
        iter <- calcobj(X, iter, pa)
        if (iter$obj > best.iter$obj)
            best.iter = iter
        best.iter_ = list(iter = best.iter, pa = pa, X = X)
    }
    return(best.iter_)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
