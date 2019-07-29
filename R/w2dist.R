w2dist = function(P,Q){
  sqrt(abs(optimalFlow::distGaussianMean(P, Q) + optimalFlow::distGaussianCov(P, Q)))
}
