wasserCostFunction <- function(j, i, cytometries, equal.weights = FALSE){
  cyto_a = cytometries[[j]]
  cyto_b = cytometries[[i]]
  names_a = unlist(lapply(cyto_a, function(x) x$type))
  names_b = unlist(lapply(cyto_b, function(x) x$type))
  weights_a = unlist(lapply(cyto_a, function(x) x$weight))
  if(abs(1-sum(weights_a)) >10^(-5)){
    message("Weights do not add to 1. A normalization will be applied.")
    weights_a = weights_a/sum(weights_a)
  }
  weights_b = unlist(lapply(cyto_b, function(x) x$weight))
  if(abs(1-sum(weights_b)) >10^(-5)){
    message("Weights do not add to 1. A normalization will be applied.")
    weights_b = weights_a/sum(weights_a)
  }

  if(equal.weights){
    a = matrix(1/length(names_a), nrow = 1, ncol = length(names_a))
    b = matrix(1/length(names_b), nrow = 1, ncol = length(names_b))
    colnames(a) = names_a
    colnames(b) = names_b
  } else{
    a = matrix(weights_a, nrow = 1, ncol = length(names_a))
    b = matrix(weights_b, nrow = 1, ncol = length(names_b))
    colnames(a) = names_a
    colnames(b) = names_b
  }

  nn = length(a)
  mm = length(b)
  cost_ab = array(dim = c(nn, mm))

  d_w_naive = 0
  for(k in 1:nn){
    for(l in 1:mm){
      cost_ab[k, l] = optimalFlow::w2dist(cyto_a[[k]], cyto_b[[l]])
      d_w_naive = d_w_naive + cost_ab[k, l]*a[k]*b[l]
    }
  }
  t_plan = transport::transport(a, b, cost_ab)
  t_cost = 0
  for(k in 1:dim(t_plan)[1]){
    t_cost = t_cost + cost_ab[t_plan[k, 1], t_plan[k, 2]]*t_plan[k, 3]
  }
  transport_cost = t_cost/d_w_naive
}
