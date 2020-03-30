#COMPATIBILITY WITH TCLUST CRAN names of input and output parameters and their corresponding equivalences

#input parameters in TCLUST CRAN and equivalence here in brackets
#x A matrix or data (x)
#k The number of clusters (k)
#alpha The proportion of observations to be trimmed (alpha)
#nstart The number of random initializations to be performed (nstart)
#iter.max The maximum number of concentration steps to be performed (iter.max)
#restr The type of restriction to be applied on the cluster scatter matrices "eigen"  "deter" and "sigma" (restr)
#restr.fact = 1.
#equal.weights specifying whether equal cluster weights are equal (equal.weights)
#zero.tol The zero tolerance used. By default set to 1e-16 (zero.tol)


#output parameters in TCLUST CRAN and equivalence here in brackets
#centers A matrix of size p x k containing the centers columwise of each cluster. (centers)
#cov An array of size p x p x k containing the covariance matrices of each cluster (cov)
#cluster A numerical vector of size n containing the cluster assignment for each observation (clusters)
#obj The value of the objective function of the best solution (obj)
#weights A numerical vector of length k, containing the weights of each cluster (weights)


#####FUNCTIONS DEVOTED TO APPLY CONSTRAINTS TO COVARIANCE MATRICES : restr2_eigen  restr2_deter_ restr.avgcov restr.diffax

##restr2_eigen
##FUNCTION FOR APPLYING EIGEN CONSTRAINTS. These are the typical constraints
## Fritz, H., Garcia-Escudero, L. A., & Mayo-Iscar, A. (2012). tclust: An r package for a trimming approach to cluster analysis. Journal of Statistical Software, 47(12), 1-26.

#' restr2_eigenv
#'
#' Function for applying eigenvalue constraints. These are the typical constraints used in tclust.
#'
#' @param autovalues Matrix containing eigenvalues.
#' @param ni.ini Current sample size of the clusters.
#' @param factor_e The level of the constraints.
#' @param zero.tol Tolerance level.
#'
#' @return The restricted eigenvalues.
#'
#' @examples
#' autovalues <- matrix(c(2, 3, 4, 1, 2, 3), nrow = 2)
#' ni.ini <- c(2, 2, 3)
#' factor_e <- 1.1
#' zero.tol <- 1e-9
#' autovalues_const <- restr2_eigenv(autovalues, ni.ini, factor_e, zero.tol)
#' autovalues_const
#'
#' @noRd
#'
restr2_eigenv <- function(autovalues, ni.ini, factor_e, zero.tol) {
###### function parameters:
###### autovalues: matrix containin eigenvalues
###### ni.ini: current sample size of the clusters
###### factor_e: the level of the constraints
###### zero.tol:toletance level

  	###### inicializations
    c <- factor_e
  	d <- t(autovalues)
    p <- nrow(autovalues)
  	K <- ncol(autovalues)
  	n <- sum(ni.ini)
  	nis <- matrix(data = ni.ini, nrow = K, ncol = p)

  	###### d_ is the ordered set of values in which the restriction objective function change the definition
  	###### points in d_ correspond to  the frontiers for the intervals in which this objective function has the same definition
  	###### ed is a set with the middle points of these intervals

  	d_ <- sort(c(d,d/c))
  	dim <- length(d_)
  	d_1 <- d_
  	d_1[dim + 1] <- d_[dim] * 2
  	d_2 <- c(0, d_)
  	ed <- (d_1 + d_2)/2
  	dim <- dim + 1;

  	###### the only relevant eigenvalues are those belong to a clusters with sample size greater than 0.
  	###### eigenvalues corresponding to a clusters whit 0 individuals has no influence in the objective function
  	###### if all the eigenvalues are 0 during the smart initialization we assign to all the eigenvalues the value 1

  	if ((max(d[nis > 0]) <= zero.tol))
  		  return(matrix(0, nrow = p, ncol = K))	        ##  solution corresponds to 0 matrix

  	###### we check if the  eigenvalues verify the restrictions

  	if (abs(max(d[nis > 0])/min(d[nis > 0])) <= c) {
    		d[nis == 0] <- mean(d[nis > 0])
    		return(t(d))					## the solution correspond to the input becuse it verifies the constraints
    		#dfin=d
  	}

  	###### we build the sol array
  	###### sol[1],sol[2],.... this array contains the critical values of the interval functions which defines the m objective function
  	###### we use the centers of the interval to get a definition for the function in each interval
  	###### this set with the critical values (in the array sol) contains the optimum m value

  	t <- s <- r <- array(0, c(K, dim))
  	sol <- sal <- array(0, c(dim))

  	for (mp_ in 1:dim) {
    		for (i in 1:K) {
      			r[i, mp_] <- sum((d[i, ] < ed[mp_])) + sum((d[i, ] > ed[mp_] * c))
      			s[i, mp_] <- sum(d[i, ] * (d[i, ] < ed[mp_]))
      			t[i, mp_] <- sum(d[i, ] * (d[i, ] > ed[mp_] * c))
    		}

    		sol[mp_] <- sum(ni.ini/n * (s[, mp_] + t[, mp_]/c))/(sum(ni.ini/n * (r[, mp_])))

    		e <- sol[mp_] * (d < sol[mp_]) + d * (d >= sol[mp_])*(d <= c * sol[mp_]) + (c * sol[mp_]) * (d > c * sol[mp_])
    		o <- -1/2 * nis/n * (log(e) + d/e)

    		sal[mp_] <- sum(o)
  	}

  	###### m is the optimum value for the eigenvalues procedure
  	eo <- which.max(c(sal))
  	m <- sol[eo]

  	###### based on the m value we get the restricted eigenvalues

  	t(m * (d < m) + d * (d >= m) * (d <= c * m) + (c * m) * (d > c * m))	##	the return value
}


#      restr.diffax
##	function which manages constraints application

#' restr.diffax
#'
#' Ment for applying restrictions to solutions obtained by tclust_.
#'
#' @param iter Current solution obtained by tclust_.
#' @param pa Parameters for using in tclust_.
#'
#' @return The new restricted solution.
#'
#' @examples
#' x <- rbind(matrix(rnorm(100), ncol = 2), matrix(rnorm(100) + 2, ncol = 2), matrix(rnorm(100) + 4, ncol = 2))
#' output3 <- tclust_(X = x, K = 3, alpha = 0.05, niter = 20, Ksteps = 10,
#'                  equal.weights = FALSE, restr.cov.value = "eigen",
#'                 maxfact_e = 5, zero.tol = 1e-16, trace = 0,
#'                 sol_ini_p = FALSE, sol_ini = NA)
#' iter <- output3$iter
#' pa <- output3$pa
#' pa$maxfact_e <- 1.1
#' print(iter$sigma)
#' iter2 <- restr.diffax(iter, pa)
#' print(iter2$sigma)
#'
#' @noRd
#'
restr.diffax <- function (iter, pa) {

  	u <- array(NA, c(pa$p, pa$p, pa$K))
  	d <- array(NA, c(pa$p, pa$K))

  	for (k in 1:pa$K) {
    		ev <- eigen(iter$sigma[, , k])
    		u[, , k] <- ev$vectors
    		d[, k] <- ev$values
  	}

  	d[d < 0] <- 0		##	all eigenvalue < 0 are assigned to 0, this issue appears for numerical errors

    d <- restr2_eigenv(autovalues = d, ni.ini = iter$csize, factor_e = pa$maxfact_e, zero.tol = pa$zero.tol)
  	##	checking for singularity in all clusters.
  	iter$code <- max(d) > pa$zero.tol

  	if (!iter$code)
    		return(iter)

  	for (k in 1:pa$K)	##	re-composing the sigmas
  		  iter$sigma[, , k] <- u[, , k] %*% diag(d[, k], nrow = pa$p) %*% t(u[, , k])

  	return(iter)
}




###### MISCELANEOUS FUNCTIONS: dmnorm ssclmat TreatSingularity

##	Multivariate normal density

#' dmnorm
#'
#' The multivariate normal density.
#'
#' @param X Multivariate points.
#' @param mu Mean.
#' @param sigma Covariance matrix.
#'
#' @return The value of the multivariate normal with mean mu and covariance sigma at the points in X.
#' @examples
#' x <- rbind(matrix(rnorm(100), ncol = 2), matrix(rnorm(100) + 2, ncol = 2), matrix(rnorm(100) + 4, ncol = 2))
#' dmnorm(X = x, mu = c(0, 0), sigma = diag(2))
#'
#' @importFrom stats mahalanobis
#'
#' @noRd
#'
dmnorm <- function (X, mu, sigma) ((2 * pi)^(-length(mu)/2)) * (det(sigma)^(-1/2)) * exp(-0.5 * mahalanobis(X, mu, sigma))


##	get a matrix object out of the sigma

#' ssclmat
#'
#' Extract a matrix from the object containing covariance matrices. Ment for use inside tclust_ function.
#'
#' @param x An object used inside tclust_ containing covariance matrices.
#' @param k An integer value giving a location in an array.
#'
#' @return A covariance matrix.
#' @examples
#' x <- rbind(matrix(rnorm(100), ncol = 2), matrix(rnorm(100) + 2, ncol = 2), matrix(rnorm(100) + 4, ncol = 2))
#' output3 <- tclust_(X = x, K = 3, alpha = 0.05, niter = 20, Ksteps = 10,
#'                  equal.weights = FALSE, restr.cov.value = "eigen",
#'                  maxfact_e = 5, zero.tol = 1e-16, trace = 0,
#'                  sol_ini_p = FALSE, sol_ini = NA)
#' iter <- output3$iter
#' pa <- output3$pa
#' pa$maxfact_e <- 1.1
#' iter2 <- restr.diffax(iter, pa)
#' ssclmat(iter2$sigma, k = 1)
#'
#' @noRd
#'
ssclmat <- function (x, k) as.matrix(x[, , k])

##	to manage singular situations

#' TreatSingularity
#'
#' Shows a warning message when a singularity in the curent solution of tclust_ is found.
#'
#' @param iter Current solution obtained by tclust_.
#' @param pa Parameters for using in tclust_.
#'
#' @return A warning message.
#' @examples
#' x <- rbind(matrix(rnorm(100), ncol = 2), matrix(rnorm(100) + 2, ncol = 2), matrix(rnorm(100) + 4, ncol = 2))
#' output3 <- tclust_(X = x , K = 3, alpha = 0.05, niter = 20, Ksteps = 10,
#'                  equal.weights = FALSE, restr.cov.value = "eigen",
#'                  maxfact_e = 5, zero.tol = 1e-16, trace = 0,
#'                  sol_ini_p = FALSE, sol_ini = NA)
#' iter <- output3$iter
#' pa <- output3$pa
#' TreatSingularity(iter, pa)
#'
#' @noRd
#'
TreatSingularity <- function (iter, pa) {
  	warning("points in the data set are concentrated in k points after trimming ")
  	return(iter)
}


######## FUNCTIONS FOR RANDOM STARTS:  getini    InitClusters
##	calculates the initial cluster sizes

#' getini
#'
#' Calculates the initial cluster sizes for initializing tclust_.
#'
#' @param K Number of groups.
#' @param no.trim Indicates if trimming is allowed.
#'
#' @return A vector of length K indicating indices of the points to be considered for initialization.
#' @examples
#' v <- getini(K = 3, no.trim = 100)
#' v
#'
#' @noRd
#'
getini <- function (K, no.trim) {
  	if (K == 1)
  		  return(no.trim)

  	pi.ini <- runif(K)
  	ni.ini <- sample(x = K, size = no.trim, replace = TRUE, prob = pi.ini/sum(pi.ini))
  	return(tabulate(ni.ini, nbins = K))
}

##	calculates the initial cluster assignment and initial values for the parameters

#' InitClusters
#'
#' Calculates the initial cluster assignment and initial values for the parameters
#'
#' @param X Points in a multidimiensional space.
#' @param iter Current solution obtained by tclust_.
#' @param Parameters for using in tclust_.
#'
#' @return An initial solution, based on a random subsample, in a form of an object used internally in tclust_.
#' @examples
#' x <- rbind(matrix(rnorm(100), ncol = 2), matrix(rnorm(100) + 2, ncol = 2), matrix(rnorm(100)+ 4, ncol = 2))
#' output3 <- tclust_(X = x, K = 3, alpha = 0.05, niter = 20, Ksteps = 10,
#'                  equal.weights = FALSE, restr.cov.value = "eigen",
#'                  maxfact_e = 5, zero.tol = 1e-16,  trace = 0,
#'                  sol_ini_p = FALSE, sol_ini = NA)
#' iter <- output3$iter
#' pa <- output3$pa
#' iter <- InitClusters(X = x, iter = output3$iter, pa = output3$pa)
#' iter$cw
#' iter$center
#' iter$sigma
#'
#' @importFrom stats cov
#'
#' @noRd
#'
InitClusters <- function (X, iter, pa) {
  	dMaxVar <- 0
  	for (k in 1:pa$K) {
    		idx <- sample(1:pa$n, pa$p + 1)
    		X.ini <- X[drop = FALSE, idx, ]#sample (1:pa$n, pa$p+1),]	##	selecting observations randomly for the current init - cluster

    		iter$center[k, ] <- colMeans(X.ini)			##	calculating the center

    		cc <- (pa$p/(pa$p + 1)) * cov(X.ini)			##	calculating sigma (cc = current cov)
    		iter$sigma[, , k] <- cc
  	}
  	if (pa$equal.weights) {						##	if we're considering equal weights, cw is set here AND NEVER CHANGED
  	    iter$csize <- rep(pa$no.trim/pa$K, pa$K)
  	}	else {
  	      iter$csize = getini(pa$K, pa$no.trim)
  	}
  	iter$cw <- iter$csize/pa$no.trim				## if we're considering different weights, calculate them, and they're gonna be recalculated every time in findClustAss
  	return(iter)
}


######## FUNCTION FOR estimating model parameters: estimClustPar
##PAPERS
##HARD ASSIGNMEMT Fritz, H., Garcia-Escudero, L. A., & Mayo-Iscar, A. (2012). tclust: An r package for a trimming approach to cluster analysis. Journal of Statistical Software, 47(12), 1-26.
##MIXTURE MODEL Garc?a-Escudero, Luis Angel, Alfonso Gordaliza, and Agust?n Mayo-Iscar. "A constrained robust proposal for mixture modeling avoiding spurious solutions." Advances in Data Analysis and Classification 8.1 (2014): 27-43.
##FUZZY CLUSTERING Fritz, H., Garc?A-Escudero, L. A., & Mayo-Iscar, A. (2013). Robust constrained fuzzy clustering. Information Sciences, 245, 38-52.

#' estimClustPar
#'
#' Obtain the best values for the model parameters, given data, input parameters and an assigment.
#'
#' @param X Points in a multidimiensional space.
#' @param iter Current solution obtained by tclust_.
#' @param pa Parameters for using in tclust_.
#'
#' @return An object used internally in tclust_ with the best model parameters.
#' @examples
#' x <- rbind(matrix(rnorm(100), ncol = 2), matrix(rnorm(100) + 2, ncol = 2), matrix(rnorm(100) + 4, ncol = 2))
#' output3 <- tclust_(X = x, K = 3, alpha = 0.05, niter = 20, Ksteps = 10,
#'                  equal.weights = FALSE, restr.cov.value = "eigen",
#'                  maxfact_e = 5, zero.tol = 1e-16, trace = 0,
#'                  sol_ini_p = FALSE, sol_ini = NA)
#' iter <- output3$iter
#' pa <- output3$pa
#' output4 <- estimClustPar(X = x, iter, pa)
#' output4$center
#' output4$sigma
#' output4$cw
#'
#' @noRd
#'
estimClustPar <- function (X, iter, pa) {
    ###  this line must dissapear because  iter$csize  is updated in findClustAssig          iter$csize=apply(iter$z_ij,2,'sum')
		for (k in 1:pa$K) {
  			if (iter$csize[k] > pa$zero.tol) { ##	if cluster's size is > 0
    				iter$center[k, ] = (t(iter$z_ij[, k]) %*% X)/iter$csize[k]
    				X.c <- (X - matrix(iter$center[k, ], ncol = pa$p, nrow = pa$n, byrow = TRUE))
    				iter$sigma[, , k] <- (t(X.c * iter$z_ij[, k]) %*% X.c)/iter$csize[k]
  			}	else							##	this cluster's size has decreased to 0
  			    iter$sigma[, , k] <- 0
    }
	return(iter)
}


######## FUNCTIONS FOR obtaining the assigment and trimming: findClustAssig (mixture models and hard assigment)  findClustAssig_  (fuzzy assigment)

######## FUNCTION FOR obtaining the assigment and trimming in the non FUZZY CASE (mixture and hard assignments)
##PAPERS
##HARD ASSIGNMEMT Fritz, H., Garcia-Escudero, L. A., & Mayo-Iscar, A. (2012). tclust: An r package for a trimming approach to cluster analysis. Journal of Statistical Software, 47(12), 1-26.
##MIXTURE MODEL Garc?a-Escudero, Luis Angel, Alfonso Gordaliza, and Agust?n Mayo-Iscar. "A constrained robust proposal for mixture modeling avoiding spurious solutions." Advances in Data Analysis and Classification 8.1 (2014): 27-43.

#' findClustAssig
#'
#' Function for obtaining the assigment and trimming (mixture models and hard assigment).
#'
#' @param X Points in a multidimiensional space.
#' @param iter Current solution obtained by tclust_.
#' @param pa Parameters for using in tclust_.
#'
#' @return An object used internally in tclust_ with the assignment of the input points.
#' @examples
#' x <- rbind(matrix(rnorm(100), ncol = 2), matrix(rnorm(100) + 2, ncol = 2), matrix(rnorm(100) + 4, ncol = 2))
#' output3 <- tclust_(X = x, K = 3, alpha = 0.05, niter = 20, Ksteps = 10,
#'                  equal.weights = FALSE, restr.cov.value = "eigen",
#'                  maxfact_e = 5, zero.tol = 1e-16, trace = 0,
#'                  sol_ini_p = FALSE, sol_ini = NA)
#' iter <- output3$iter
#' pa <- output3$pa
#' output5 <- findClustAssig(X=x, iter, pa)
#'
#' @noRd
#'
findClustAssig <- function (X, iter, pa) {
    ll = matrix(NA, pa$n, pa$K)
    for (k in 1:pa$K)	ll[, k] <- iter$cw[k]*dmnorm(X, iter$center[k, ], ssclmat(iter$sigma, k)) ## dmvnorm could be used here...

  	old.assig <- iter$assig
  	iter$assig <- apply(ll, 1, which.max)						##	searching the cluster which fits best for each observation

    pre.z_h <- apply(ll, 1, 'max')

    pre.z_m <- apply(ll, 1, 'sum')

    pre.z_ <- matrix(pre.z_m, nrow = pa$n, ncol = pa$K, byrow = FALSE)

    ##### To obtain the trimming:  tc.set is the non trimming indicator

    tc.set <- (rank(pre.z_h, ties.method = "random") > floor(pa$n * (pa$alpha)))

    ##### To obtain the iter$z_ij matrix contining the assigment and trimming
    iter$assig <- apply(ll, 1, which.max)*tc.set    #hard assigment including trimming

    iter$z_ij <- 0 * iter$z_ij
    iter$z_ij[cbind((1:pa$n), iter$assig + (iter$assig == 0))] <- 1
    iter$z_ij[tc.set == FALSE, ] <- 0
    iter$code <- 2 * all(old.assig == iter$assig)		##	setting the code - parameter, signaling whether the assignment is the same than the previous --- is the only stopping rule implemented

    ##### To obtain the size of the clusters and the estimated weight of each population
  	iter$csize <- tabulate(iter$assig, pa$K)

    if (!pa$equal.weights) iter$cw <- iter$csize/sum(iter$csize)					##	and cluster weights

    return(iter)
}




######## FUNCTION FOR obtaining the objetive functions value for mixture (obj_m) hard (obj_h) and fuzzy (obj_f)
##PAPERS
##HARD ASSIGNMEMT Fritz, H., Garcia-Escudero, L. A., & Mayo-Iscar, A. (2012). tclust: An r package for a trimming approach to cluster analysis. Journal of Statistical Software, 47(12), 1-26.
##MIXTURE MODEL Garc?a-Escudero, Luis Angel, Alfonso Gordaliza, and Agust?n Mayo-Iscar. "A constrained robust proposal for mixture modeling avoiding spurious solutions." Advances in Data Analysis and Classification 8.1 (2014): 27-43.
##FUZZY CLUSTERING Fritz, H., Garc?A-Escudero, L. A., & Mayo-Iscar, A. (2013). Robust constrained fuzzy clustering. Information Sciences, 245, 38-52.

#' calcobj
#'
#' Calculates tclust's objective function value.
#'
#' @param X Points in a multidimiensional space.
#' @param iter Current solution obtained by tclust_.
#' @param pa Parameters for using in tclust_.
#'
#' @return An object used internally in tclust_ with the the value for tclust's objective function.
#' @examples
#' x <- rbind(matrix(rnorm(100), ncol = 2), matrix(rnorm(100) + 2, ncol = 2),
#'         matrix(rnorm(100) + 4, ncol = 2))
#' output3 <- tclust_(X = x, K = 3, alpha = 0.05, niter = 20, Ksteps = 10,
#'                  equal.weights = FALSE, restr.cov.value = "eigen",
#'                  maxfact_e = 5, zero.tol = 1e-16, trace = 0,
#'                  sol_ini_p = FALSE, sol_ini = NA)
#' iter <- output3$iter
#' pa <- output3$pa
#' iter_ <- calcobj(X = x, iter = iter, pa = pa)
#' iter_$obj
#'
#' @noRd
#'
calcobj <- function (X, iter, pa) {
    ww_m <- matrix(0, nrow = pa$n, ncol = 1)
    ww_h <- matrix(0, nrow = pa$n, ncol = 1)

  	for (k in 1:pa$K) {
         w_m <- iter$cw[k] * dmnorm(X, iter$center[k, ], ssclmat(iter$sigma, k))
         ww_m <- w_m * (w_m >= 0) + ww_m     ##	calculates each individual contribution for the obj funct mixture
         w_h <- w_m * (iter$assig == k)
         ww_h <- w_h * (w_h >= 0) + ww_h     ##	calculates each individual contribution for the obj funct hard
    }
    ww_m <- ww_m * (ww_m >= 0)
    ww_h <- ww_h * (ww_h >= 0)

    iter$obj <- sum(log(ww_h[iter$assig > 0]))
  	return(iter)
}
#' tclust_
#'
#' Function that performs robust non spherical clustering, tclust, where initial solutions are allowed.
#'
#' @param X Points in a multidimiensional space.
#' @param K Number of clusters.
#' @param alpha Level of trimming.
#' @param niter Maximum number of iterations.
#' @param Ksteps Maximum number of K-steps.
#' @param equal.weights If all clusters should have equal weight in the objective function.
#' @param restr.cov.value Type of restriction on covariance matrices. Only "eigen" will be allowed.
#' @param maxfact_e Level determinant constraints.
#' @param zero.tol The zero tolerance used. By default set to 1e-16.
#' @param trace Defines the tracing level, which is set to 0 by default. Tracing level 2 gives additional information on the iteratively decreasing objective function's value.
#' @param sol_ini_p Initial solution for parameters provided by the user TRUE/FALSE, if TRUE is stored in sol_ini.
#' @param sol_ini Initial solution for parameters provided by the user.
#'
#' @return A list with values:
#' \describe{
#'  \item{iter}{Contains the solutions of the clustering provided by tclust_.}
#'  \item{pa}{Parametaers used in tslucs_.}
#'  \item{X}{The input variables.}
#' }
#'
#' @examples
#' x <- rbind(matrix(rnorm(100), ncol = 2), matrix(rnorm(100) + 2, ncol = 2), matrix(rnorm(100) + 4, ncol = 2))
#' output3 <- tclust_(X = x, K = 3, alpha = 0.05, niter = 20, Ksteps = 10,
#'                  equal.weights = FALSE, restr.cov.value = "eigen",
#'                  maxfact_e = 5, zero.tol = 1e-16, trace = 0,
#'                  sol_ini_p = FALSE, sol_ini = NA)
#' output3$iter$assig
#' plot(x, col = output3$iter$assig + 1)
#'
#' @noRd
#'
tclust_ <- function(X, K, alpha = 0.05, niter = 20, Ksteps = 10, equal.weights = FALSE,
                    restr.cov.value = "eigen", maxfact_e = 5, zero.tol = 1e-16, trace = 0,
                    sol_ini_p = FALSE, sol_ini = NA) {

	  if (!is.numeric(X))
		    stop("parameter x: numeric matrix/vector expected")
	  if (!is.matrix(X))
		    X <- matrix(X, ncol = 1)
  	n <- nrow(X)
  	p <- ncol(X)
  	no.trim <- floor(n * (1 - alpha))
    f.restr <- restr.diffax

  	# preparing lot's of lists: pa ( input parameters for the procedure)  iter (current value of the parameters of the model)
  	pa <- list(		            ##	 input parameters for the procedure
    		n = n,			    ##	number of observations
    		p = p,	                    ##	number of dimensions
    		alpha = alpha,		    ##	level of trimming
    		trimm = n - no.trim,	    ##	number of observations which are considered as to be outliers
    		no.trim = no.trim,	    ##	number of observations which are considered as to be not outliers
    		K = K,			    ##	number of clusters
    		equal.weights = equal.weights,		##	equal population proportions  for all clusters
    		zero.tol = zero.tol,			##	zero tolerance	(to be implemented)
    		trace = trace,				##	level of information provided (to be implemented)
        maxfact_e = maxfact_e,                    ##	level determinant constraints
        sol_ini_p = sol_ini_p,                  ##	initial solution for parameters provided by the user TRUE/FALSE   if TRUE is stored in sol_ini
        sol_ini = sol_ini,                      ##	initial solution for parameters provided by the user
        Ksteps = Ksteps,                          ##	number of iterations
        niter = niter                             ##	number of random starts
  	)

  	iter <- list(					##	current value of the parameters
    		obj = NA,				##	current objective value
    		assig = array(0, n),			##	cluster group of assignment
    		csize = array(NA, K),			##	cluster number of observations assigned to the cluster
    		cw = rep(NA, K),			##	cluster estimated proportions
    		sigma = array(NA, c (p, p, K)),	##	cluster's covariance matrix
    		center = array(NA, c(K, p)),		##	cluster's centers
    		code = NA,				##	this is a return code supplied by functions like findClustAssig
    		z_ij = matrix(0, nrow = n, ncol = K),	##	cluster assignment given by 0/1 columns
    		lambda = rep(NA, K)	##	diagonal values for tkmeans
  	)
	  best.iter <- list(obj = -Inf)			##	for containing the best values for the paramteters after several random starts

    if (pa$sol_ini_p == TRUE) pa$niter <- 1

	  for (j in 1:pa$niter) {
      if (pa$sol_ini_p == TRUE) {
        if (is.numeric(sol_ini$assig)) {
            iter$assig <- sol_ini$assig
            iter$z_ij <- 0 * iter$z_ij
            for (h in 1:pa$K) iter$z_ij[iter$assig == h, h] <- 1
  	        iter$csize <- tabulate(iter$assig, pa$K)
            iter$cw <- iter$csize/sum(iter$csize)
   	        iter <- estimClustPar(X, iter, pa)
   	        iter$center[is.na(iter$center)] <- 0
        } else  {
            iter$cw <- pa$sol_ini$cw;  iter$center <- pa$sol_ini$center; iter$sigma <- pa$sol_ini$sigma
	          iter$csize <- iter$cw*pa$n * (1 - pa$alpha)
	      }
      } else  iter <- InitClusters(X, iter, pa)   #initial solution provided by the user (pa$sol_ini_p==TRUE) or at random (FALSE)

      iter$center[is.na(iter$center)] <- 0
		  for (i in 0:pa$Ksteps) {
          iter <- f.restr(iter = iter, pa = pa)	##	restricting the clusters' scatter structure
    			if (iter$code == 0) {								##	all eigenvalues are zero
    				if (i > 0)
    					return(TreatSingularity(calcobj(X, iter, pa), pa))
    				else
    					iter$sigma[, , ] <- diag(pa$p)
          }
    		  iter <- findClustAssig(X, iter, pa)		## estimates the cluster's assigment and TRIMMING (mixture models and HARD )
    			if ((iter$code == 2) ||									##		if findClustAssig returned 1, meaning that the cluster Assignment has not changed
    		         (i == pa$Ksteps))									##		or we're in the last concentration step:
    				break											##		break the for - loop - we finished this iteration! dont re-estimate cluster parameters this time
    			iter <- estimClustPar(X, iter, pa)		## estimates the cluster's parameters
		  }
		  iter <- calcobj(X, iter, pa)							##	calculates the objetive function value
		if (iter$obj > best.iter$obj)
			best.iter <- iter
    best.iter_ <- list(iter = best.iter, pa = pa, X = X)
	}
	return(best.iter_)
}
#' tclust_H
#'
#' A wrapper for the internal fucntion tclust_. Performs robust non spherical clustering, tclust, where initial solutions are allowed.
#'
#' @param x A matrix or data.frame of dimension n x p, containing the observations (row-wise).
#' @param k The number of clusters initially searched for.
#' @param alpha The proportion of observations to be trimmed.
#' @param nstart The number of random initializations to be performed. Only when sol_ini_p = FALSE.
#' @param iter.max The maximum number of concentration steps to be performed. The concentration steps are stopped, whenever two consecutive steps lead to the same data partition.
#' @param restr The type of restriction to be applied on the cluster scatter matrices. Valid values are "eigen" (default).
#' @param restr.fact The constant restr.fact >= 1 constrains the allowed differences among group scatters. Larger values imply larger differences of group scatters, a value of 1 specifies the strongest restriction.
#' @param sol_ini_p Initial solution for parameters provided by the user TRUE/FALSE, if TRUE is stored in sol_ini.
#' @param sol_ini Initial solution for parameters provided by the user.
#' @param equal.weights A logical value, specifying whether equal cluster weights (TRUE) or not (FALSE) shall be considered in the concentration and assignment steps.
#' @param trace Defines the tracing level, which is set to 0 by default. Tracing level 2 gives additional information on the iteratively decreasing objective function's value.
#' @param zero.tol The zero tolerance used. By default set to 1e-16.
#'
#' @return A list with values:
#' \describe{
#'  \item{centers}{A matrix of size p x k containing the centers (column-wise) of each cluster.}
#'  \item{cov}{An array of size p x p x k containing the covariance matrices of each cluster.}
#'  \item{cluster}{A numerical vector of size n containing the cluster assignment for each observation. Cluster names are integer numbers from 1 to k, 0 indicates trimmed observations.}
#'  \item{par}{A list, containing the parameters the algorithm has been called with (x, if not suppressed by store.x = FALSE, k, alpha, restr.fact, nstart, KStep, and equal.weights).}
#'  \item{weights}{A numerical vector of length k, containing the weights of each cluster.}
#'  \item{obj}{he value of the objective function of the best (returned) solution.}
#' }
#' @details
#' This iterative algorithm initializes k clusters randomly and performs "concentration steps" in order to improve the current cluster assignment. The number of maximum concentration steps to be performed is given by iter.max. For approximately obtaining the global optimum, the system is initialized nstart times and concentration steps are performed until convergence or iter.max is reached. When processing more complex data sets higher values of nstart and iter.max have to be specified (obviously implying extra computation time). However, if more then half of the iterations would not converge, a warning message is issued, indicating that nstart has to be increased.
#'
#' The parameter restr defines the cluster's shape restrictions, which are applied on all clusters during each iteration. Options "eigen"/"deter" restrict the ratio between the maximum and minimum eigenvalue/determinant of all cluster's covariance structures to parameter restr.fact. Setting restr.fact to 1, yields the strongest restriction, forcing all eigenvalues/determinants to be equal and so the method looks for similarly scattered (respectively spherical) clusters. Option "sigma" is a simpler restriction, which averages the covariance structures during each iteration (weighted by cluster sizes) in order to get similar (equal) cluster scatters.
#'
#' @examples
#' x <- rbind(matrix(rnorm(100), ncol = 2), matrix(rnorm(100) + 2, ncol = 2),
#'         matrix(rnorm(100) + 4, ncol = 2))
#' ## robust cluster obtention from a sample x asking for 3 clusters,
#' ## trimming level 0.05 and constrain level 12
#' k <- 3; alpha <- 0.05; restr.fact <- 12
#' output <- tclust_H(x = x, k = k, alpha = alpha, nstart = 50, iter.max = 20,
#'                  restr = "eigen", restr.fact = restr.fact, sol_ini_p = FALSE, sol_ini = NA,
#'                  equal.weights = FALSE, trace = 0, zero.tol = 1e-16)
#' ## cluster assigment
#' output$cluster
#' plot(x, col = output$cluster)
#'
#' @references Fritz, H., Garcia-Escudero, L. A., & Mayo-Iscar, A. (2012). tclust: An r package for a trimming approach to cluster analysis. Journal of Statistical Software, 47(12), 1-26.
#'
#' @export
#'
tclust_H <- function (x, k = 3, alpha = 0.05, nstart = 50, iter.max = 20, restr = "eigen", restr.fact = 12,
                    sol_ini_p = FALSE, sol_ini = NA, equal.weights = FALSE, trace = 0, zero.tol = 1e-16) {
  if (sol_ini_p == TRUE) {
    if (is.numeric(sol_ini$cluster)) {
      sol_ini$assig = sol_ini$cluster
      } else {
      sol_ini$cw <- sol_ini$weights; sol_ini$center <- t(sol_ini$centers); sol_ini$sigma <- sol_ini$cov
    }
  }
  tc_ <- tclust_(X = x, K = k, alpha = alpha, niter = nstart, Ksteps = iter.max, equal.weights = equal.weights, restr.cov.value = restr,
                 maxfact_e = restr.fact,  zero.tol = 1e-16,  trace = 0,
                 sol_ini_p = sol_ini_p, sol_ini = sol_ini )
  tc <- list()
  tc$centers <- t(tc_$iter$center)
  tc$cov <- tc_$iter$sigma
  tc$cluster <- tc_$iter$assig
  tc$par <- tc_$pa
  tc$weights <- tc_$iter$cw
  tc$obj <- tc_$iter$obj
  return(tc)
}
