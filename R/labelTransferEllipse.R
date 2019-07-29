labelTransferEllipse = function(i, test.cytometry.ellipses, training.cytometries.barycenter, equal.weights = FALSE){
  N = length(training.cytometries.barycenter)
  M = length(test.cytometry.ellipses)
  cost.matrix = array(dim = c(N,M))
  for (j in 1:N){
    for(i in 1:M){
      cost.matrix[j, i] = w2dist(list(mean = test.cytometry.ellipses[[i]]$mean, cov = test.cytometry.ellipses[[i]]$cov),
                                 list(mean = training.cytometries.barycenter[[j]]$mean,
                                      cov = training.cytometries.barycenter[[j]]$cov))
    }
  }
  names.a = unlist(lapply(training.cytometries.barycenter, function(x) x$type))
  names.b = unlist(lapply(test.cytometry.ellipses, function(x) x$type))

  if(equal.weights){
    A = matrix(1/N, nrow = 1, ncol = N)
    B = matrix(1/M, nrow = 1, ncol = M)
    names(A) = names.a
    names(B) = names.b
  } else{
    A = unlist(lapply(training.cytometries.barycenter, function(x) x$weight))
    A = A/sum(A)
    names(A) = names.a
    B = unlist(lapply(test.cytometry.ellipses, function(x) x$weight))
    names(B) = names.b
  }
  optimal.transport.form.A.to.B = transport::transport(a = A, b = B, costm = cost.matrix)
  optimal.transport.form.A.to.B$to = names(B)[optimal.transport.form.A.to.B$to]
  optimal.transport.form.A.to.B$from = names(A)[optimal.transport.form.A.to.B$from]
  return(optimal.transport.form.A.to.B)
}
