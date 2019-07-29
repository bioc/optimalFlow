tclustWithInitialization = function(initialization, cytometry, i.sol.type = "points", trimming = 0.05, restr.fact = 1000){
  t0 = Sys.time()
  dim.cyto = dim(cytometry)[2]
  if (i.sol.type == "points"){
    cosa = data.frame(initialization[, 1:dim.cyto], cluster = dplyr::pull(initialization, (dim.cyto + 1)))
    n_clus = length(levels(cosa$cluster))
    levels(cosa$cluster) = 1:n_clus

    cosa_no_noise = cosa[which(as.double(cosa$cluster) > 0),]

    sol_inicial_0 = lapply(1:n_clus, estimCovCellGeneral, cosa[, 1:dim.cyto], cosa$cluster)
    not_valid_estim = which(is.na(sol_inicial_0))
    if(length(not_valid_estim) > 0){
      sol_inicial_0 = sol_inicial_0[-not_valid_estim]
    } else{
      sol_inicial_0 = sol_inicial_0
    }
    weight_citos = unlist(lapply(sol_inicial_0, function(x) x$weight))
    means_citos = matrix(unlist(lapply(sol_inicial_0, function(x) x$mean)), ncol = length(weight_citos), byrow = FALSE)
    # cov_citos = lapply(sol_inicial_0, function(x) x$cov)
    cov_citos = array(dim = c(dim.cyto, dim.cyto, length(weight_citos)))
    for( i in 1:length(weight_citos)){
      cov_citos[,,i] = sol_inicial_0[[i]]$cov
    }

    sol_inicial = list()
    sol_inicial$weights = weight_citos
    sol_inicial$cov = cov_citos
    sol_inicial$centers = means_citos
    # sol_inicial$weights = as.vector(table(cosa_no_noise$cluster)/sum(table(cosa_no_noise$cluster)))
    KK = length(sol_inicial$weights)
    print(paste("tclust looking for k = ", KK, sep = ""))
    # covars_inicial = array(dim = c(10, 10, KK))
    # means_inicial = array(dim = c(10, KK))
    # # det_inicial = array(dim = c(10, KK))
    # for (i in 1:KK){
    #   d_cov = (cosa_no_noise %>% dplyr::filter(cluster == i))[,-11]
    #   covars_inicial[,, i] = cov(d_cov)
    #   means_inicial[, i] = colMeans(d_cov)
    #   # det_inicial[,i] = eigen(covars_inicial[,,i])$values
    # }
    # sol_inicial$cov = covars_inicial
    # sol_inicial$centers = means_inicial
    trimming = 0.5*exp(-(KK - 1)/5)
  } else{
    if(i.sol.type == "barycenter"){
      sol_inicial = initialization
      KK = length(sol_inicial$weights)
      print(paste("tclust looking for k = ", KK, sep = ""))
    }
  }

  # save(sol_inicial, file = paste(paste("sol_inicial", paste(iii, kkk, sep = "_"), sep = ""), "RData", sep = "."))
  t00 = Sys.time()
  solution.tclust = optimalFlow::tclust_H(as.matrix(cytometry[, 1:dim.cyto]), k = KK, alpha = trimming, restr.fact = restr.fact, nstart = 10, iter.max = 100, equal.weights = FALSE,
                    sol_ini = sol_inicial, sol_ini_p = TRUE, restr = "eigen",   center=FALSE, scale=FALSE, store.x = TRUE, drop.empty.clust = TRUE, trace = 0, warnings = 3,
                    zero.tol = 1e-16)
  t11 = Sys.time()
  t11-t00

  t1 = Sys.time()
  n_clus = length(table(solution.tclust$cluster)) - 1
  print(paste("tclust found k = ", n_clus, sep = ""))
  print(paste(difftime(t1, t0), units(difftime(t1, t0))))
  cat("\n")
  return(list(cluster = solution.tclust$cluster, n_clus = n_clus, obj = solution.tclust$obj))
}
