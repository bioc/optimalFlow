\name{calcobj}
\alias{calcobj}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{calcobj
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
calcobj(X, iter, pa)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{iter}{
%%     ~~Describe \code{iter} here~~
}
  \item{pa}{
%%     ~~Describe \code{pa} here~~
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
function (X, iter, pa)
{
    ww_m = matrix(0, nrow = pa$n, ncol = 1)
    ww_h = matrix(0, nrow = pa$n, ncol = 1)
    ww_f = matrix(0, nrow = pa$n, ncol = 1)
    for (k in 1:pa$K) {
        if (pa$tk == FALSE)
            w_m = iter$cw[k] * dmnorm(X, iter$center[k, ], ssclmat(iter$sigma,
                k))
        if (pa$tk == TRUE)
            w_m = iter$cw[k] * dmnorm_tk(X, iter$center[k, ],
                iter$lambda[k])
        ww_m = w_m * (w_m >= 0) + ww_m
        w_h = w_m * (iter$assig == k)
        ww_h = w_h * (w_h >= 0) + ww_h
        w_f = iter$z_ij[, k] * log(w_m * (w_m >= 0) + 1 * (iter$z_ij[,
            k] == 0))
        ww_f = w_f + ww_f
    }
    ww_m = ww_m * (ww_m >= 0)
    ww_h = ww_h * (ww_h >= 0)
    if (pa$opt == "MIXTURE")
        iter$obj <- sum(log(ww_m[iter$assig > 0]))
    if (pa$opt == "HARD")
        iter$obj <- sum(log(ww_h[iter$assig > 0]))
    if (pa$opt == "FUZZY")
        iter$obj <- sum(ww_f[iter$assig > 0])
    return(iter)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
