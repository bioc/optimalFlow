costWasserMatchingEllipse = function(test.cytometry, training.cytometries, equal.weights = FALSE){
  dim.test.cytometries = length(training.cytometries)
  cost.vector = array(0, dim = c(1, dim.test.cytometries))
  training.cytometries_elipses = list()
  for(i in 1:(dim.test.cytometries)){
    training.cytometry = training.cytometries[[i]]
    names.a = unlist(lapply(test.cytometry, function(x) x$type))
    names.b = unlist(lapply(training.cytometry, function(x) x$type))
    weights.a = unlist(lapply(test.cytometry, function(x) x$weight))
    if(abs(1-sum(weights.a)) >10^(-5)){
      message("Pesos no normalizados. Procediendo a la nor,alizacion de los pesos.")
      weights.a = weights.a/sum(weights.a)
    }
    weights.b = unlist(lapply(training.cytometry, function(x) x$weight))
    if(abs(1-sum(weights.b)) >10^(-5)){
      message("Pesos no normalizados. Procediendo a la nor,alizacion de los pesos.")
      weights.b = weights.a/sum(weights.a)
    }

    if(equal.weights){
      a = matrix(1/length(names.a), nrow = 1, ncol = length(names.a))
      b = matrix(1/length(names.b), nrow = 1, ncol = length(names.b))
      colnames(a) = names.a
      colnames(b) = names.b
    } else{
      a = matrix(weights.a, nrow = 1, ncol = length(names.a))
      b = matrix(weights.b, nrow = 1, ncol = length(names.b))
      colnames(a) = names.a
      colnames(b) = names.b
    }

    nn = length(a)
    mm = length(b)
    cost.ab = array(dim = c(nn, mm))

    naive.transport.cost = 0
    for(k in 1:nn){
      for(l in 1:mm){
        cost.ab[k, l] = optimalFlow::w2dist(test.cytometry[[k]], training.cytometry[[l]])
        naive.transport.cost = naive.transport.cost + cost.ab[k, l]*a[k]*b[l]
      }
    }
    tranport.plan = transport::transport(a, b, cost.ab)
    transport.cost = 0
    for(k in 1:dim(tranport.plan)[1]){
      transport.cost = transport.cost + cost.ab[tranport.plan[k, 1], tranport.plan[k, 2]]*tranport.plan[k, 3]
    }
    cost.vector[1,i] = transport.cost/naive.transport.cost
  }
  # cost.vector[lower.tri(cost.vector, diag = F)] = cost.vector[upper.tri(cost.vector, diag = F)]
  return(cost.vector)
}
