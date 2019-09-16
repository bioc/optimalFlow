wassersteinKBarycenter = function(i = 1, k, alpha = 0, initialization = "rnd", pooled.clusters){
  # t0 = Sys.time()

  bar = tryCatch(optimalFlow::trimmedKBarycenter(k = k, alpha0 = alpha, type.ini = initialization, pooled.clusters), error = function(x) list(bar = "error", variacion_wasser = Inf))
  variation = bar$variacion_wasser

  # t1 = Sys.time()
  # print(t1-t0)
  return(list(wasserstein.var = variation, wasserstein.k.barycenter = bar))
}
