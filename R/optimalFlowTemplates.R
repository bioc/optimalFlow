optimalFlowTemplates <- function(database, database.names = NULL, cov.estimation = "standard", alpha.cov = 0.85,
                                 equal.weights.template = TRUE, hclust.method = "complete", templates.number = NA, minPts = 2,
                                 eps = 1, consensus.method = "pooling", barycenters.number = 37, bar.repetitions = 40,
                                 alpha.bar = 0.05, bar.ini.method = "plus-plus", consensus.minPts = 3, cl.paral = 1){
  t0 = Sys.time()
  op.syst = .Platform$OS.type
  #############################################################################################################
  ############################# Step 0: Checking that entries are well defined ################################
  #############################################################################################################

  cov.estimation = match.arg(cov.estimation, c("standard", "robust"))
  if(!is.double(alpha.cov)){
    stop("alha.cov is not a doble")
  } else{
    if(alpha.cov < 0 | alpha.cov >= 1){
      stop("alpha.cov is not in the range (0,1).")
    }
  }
  hclust.method = match.arg(hclust.method, c("complete", "single", "average", "hdbscan", "dbscan"))
  if(!is.na(templates.number)){
    if(!is.double(templates.number) & !is.integer(templates.number)){
      stop("templates.number is not well defined.")
    } else{
      if(templates.number < 1){
        stop("templates.number should be >= 1.")
      }
    }
  }
  if(hclust.method %in% c("hdbscan", "dbscan")){
    minPts = as.integer(minPts)
    eps = as.double(eps)
    if(is.na(minPts) | is.na(eps)){
      stop("minPts or eps not well defined.")
    }
  }
  consensus.method = match.arg(consensus.method, c("pooling", "k-barycenter", "hierarchical"))
  if(consensus.method == "k-barycenter"){
    barycenters.number = as.integer(barycenters.number)
    if(is.na(barycenters.number)){
      stop("barycenters.number should be a number >=1.")
    }
    bar.repetitions = as.integer(bar.repetitions)
    if(is.na(bar.repetitions)){
      stop("bar.repetitions must be an integer >= 1.")
    }
    if(!is.double(alpha.bar)){
      stop("alha.bar is not a doble")
    } else{
      if(alpha.bar < 0 | alpha.bar >= 1){
        stop("alpha.bar is not in the range [0,1).")
      }
    }
    bar.ini.method = match.arg(bar.ini.method, c("rnd", "plus-plus"))
  } else{
    if (consensus.method == "hierarchical"){
      consensus.minPts = as.integer(consensus.minPts)
      consensus.eps = as.double(consensus.eps)
      if(is.na(consensus.minPts) | is.na(consensus.eps)){
        stop("consensus.minPts or consensus.eps not well defined.")
      }
    }
  }
  op.syst = match.arg(op.syst, c("unix", "windows"))

  #############################################################################################################
  ############################### Step 1: Similarity distance based clustering ################################
  #############################################################################################################
  sys.time.step1.0 = Sys.time()
  database.elliptical = list()
  for (i in 1:length(database)){
    dim.cyto = dim(database[[i]])[2]
    database.elliptical[[i]] = lapply(names(table(database[[i]][, dim.cyto])), estimCovCellGeneral, cytometry = database[[i]][,
                                                                                                                              1:(dim.cyto - 1)], labels = database[[i]][, dim.cyto], type = cov.estimation, alpha = alpha.cov)
    database.elliptical[[i]] = database.elliptical[[i]][!is.na(database.elliptical[[i]])]
  }

  # cl <- parallel::makePSOCKcluster(cl.paral)
  # doParallel::registerDoParallel(cl)
  # database.elliptical = foreach::foreach(i = 1:length(database)) %dopar% {
  #   source("estimCovCellGeneral.R")
  #   database.elliptical0 = lapply(names(table(database[[i]][, dim.cyto])), estimCovCellGeneral, cytometry = database[[i]][,
  #                              1:(dim.cyto - 1)], labels = database[[i]][, dim.cyto], type = cov.estimation, alpha = alpha.cov)
  #   database.elliptical0 = database.elliptical0[!is.na(database.elliptical0)]
  # }
  # parallel::stopCluster(cl)

  dim.cytos = length(database.elliptical)
  wasser.cost = array(0, dim = c(dim.cytos, dim.cytos))
  cytometries_elipses = list()
  if((dim.cytos - 2)*(dim.cytos - 2 + 1)/2 < cl.paral){
    n.cores = (dim.cytos - 2)*(dim.cytos - 2 + 1)/2
  } else{
    n.cores = cl.paral
  }
  cl <- parallel::makePSOCKcluster(n.cores)
  doParallel::registerDoParallel(cl)
  transport_costs_list = foreach::foreach(j = 1:(dim.cytos - 1)) %:% foreach::foreach(i = (j + 1):dim.cytos) %dopar%{
    optimalFlow::wasserCostFunction(j, i, database.elliptical, equal.weights.template)
  }
  wasser.cost[lower.tri(wasser.cost, diag = FALSE)] = unlist(transport_costs_list)
  parallel::stopCluster(cl)


  if (hclust.method %in% c("complete", "single", "average")){
    cytos.hcl = hclust(as.dist(wasser.cost), method = hclust.method)
    if(is.na(templates.number)){
      if(is.null(database.names)){
        graphics::plot(cytos.hcl, main = NULL, ylab = NULL, xlab = "")
      } else{
        graphics::plot(cytos.hcl, labels = database.names, main = NULL, ylab = NULL, xlab = "")
      }
      templates.number = as.integer(readline(prompt = "How many clusters should I look for : "))
      cytos.cluster = cutree(cytos.hcl, k = templates.number)
    } else{
      cytos.cluster = cutree(cytos.hcl, k = templates.number)
    }
  } else{
    if (hclust.method == "hdbscan"){
      cytos.cluster = dbscan::hdbscan(as.dist(wasser.cost), minPts = minPts)$cluster
    } else{
      cytos.cluster = dbscan::dbscan(as.dist(wasser.cost), minPts = minPts, eps = eps)$cluster
    }
  }
  sys.time.step1.1 = Sys.time()
  time.dif.1 = difftime(sys.time.step1.1, sys.time.step1.0)
  print(paste("step 1:", time.dif.1, units(time.dif.1), sep = " "))
  #############################################################################################################
  ####################################### Step 2: Consensus clustering ########################################
  #############################################################################################################
  sys.time.step2.0 = Sys.time()
  templates = list()
  templates.cell.types = list()
  if(consensus.method == "pooling"){ ########## pooling clusters with the same labels and taking barycenter #########
    j = 0
    for (i in sort(unique(cytos.cluster))){
      if(i == 0){
        next
      } else{
        j = j + 1
        cell.types = lapply(database.elliptical[cytos.cluster == i], function(x) lapply(x, function(y) y$type))
        cell.types = sort(unique(unlist(cell.types)))
        templates.cell.types[[j]] = cell.types
        if(length(cell.types) < cl.paral){
          n.cores = length(cell.types)
        } else{
          n.cores = cl.paral
        }
        if (op.syst == "unix"){
          templates[[j]] = parallel::mclapply(cell.types, optimalFlow::estimationCellBarycenter, database.elliptical[cytos.cluster == i],
                                              mc.cores = n.cores)
        } else{
          cl <- parallel::makePSOCKcluster(n.cores)
          templates[[j]] = parallel::parLapply(cl, cell.types, optimalFlow::estimationCellBarycenter, database.elliptical[cytos.cluster == i])
          parallel::stopCluster(cl)
        }
      }
    }
  } else{ ######## pooling all clusters for the partitions in the same group #######
    pooled.clustered.cytometries = list()
    for(i in sort(unique(cytos.cluster))){
      pooled.clustered.cytometries.0 = list()
      l = 0
      for(training.cytometry in database.elliptical[cytos.cluster == i]){
        for(k in 1:length(training.cytometry)){
          l = l + 1
          pooled.clustered.cytometries.0[[l]] = training.cytometry[[k]]
        }
      }
      pooled.clustered.cytometries[[as.character(i)]] = pooled.clustered.cytometries.0
    }
    if(consensus.method == "hierarchical"){ ######## Density based hierarchical clustering based on the wasserstein distance between
      ######## clusters. For each resulting grouping we take the one-barycenter.
      j = 0
      for (i in sort(unique(cytos.cluster))){
        if(i == 0){
          next
        } else{
          j = j + 1
          # templates.cell.types[[j]] = cell.types
          if(sum(cytos.cluster == i) == 1){
            cell.types = lapply(database.elliptical[cytos.cluster == i], function(x) lapply(x, function(y) y$type))
            cell.types = sort(unique(unlist(cell.types)))
            template.partition = cell.types
          } else{
            dim.cytos = length(pooled.clustered.cytometries[[as.character(i)]])
            wasser.cost.list = wasser.cost = array(0, dim = c(dim.cytos, dim.cytos))
            if((dim.cytos - 2)*(dim.cytos - 2 + 1)/2 < cl.paral){
              n.cores = (dim.cytos - 2)*(dim.cytos - 2 + 1)/2
            } else{
              n.cores = cl.paral
            }
            cl <- parallel::makePSOCKcluster(n.cores)
            doParallel::registerDoParallel(cl)
            wasser.cost.list = foreach::foreach(k = 1:(dim.cytos - 1)) %:% foreach::foreach(l = (k + 1):dim.cytos) %dopar%{
              optimalFlow::w2dist(pooled.clustered.cytometries[[as.character(i)]][[k]], pooled.clustered.cytometries[[
                as.character(i)]][[l]])
            }
            wasser.cost[lower.tri(wasser.cost, diag = FALSE)] = unlist(wasser.cost.list)
            parallel::stopCluster(cl)
            template.partition = dbscan::hdbscan(as.dist(wasser.cost), minPts = consensus.minPts)$cluster
          }
          template.partition.labels = sort(unique(template.partition))
          if( length(which(template.partition.labels == 0)) > 0){
            template.partition.labels = template.partition.labels[-1]
          }
          cl <- parallel::makePSOCKcluster(cl.paral)
          doParallel::registerDoParallel(cl)
          template = foreach::foreach(k = template.partition.labels) %dopar%{
            # for(k in template.partition.labels){
            pooled.elements = pooled.clustered.cytometries[[as.character(i)]][template.partition == k]
            pooled.elements.type.proportion = unlist(lapply(pooled.elements, function(x) x$type))
            pooled.elements.type.proportion = table(pooled.elements.type.proportion)
            pooled.elements.type.proportion = sort(pooled.elements.type.proportion/ sum(pooled.elements.type.proportion), decreasing = TRUE)
            pooled.elements.covs = lapply(pooled.elements, function(x) x$cov)
            pooled.elements.weighted.means = matrix(unlist(lapply(pooled.elements, function(x) x$mean)),
                                                    ncol = length(pooled.elements), byrow = FALSE)
            pooled.elements.weights = unlist(lapply(pooled.elements, function(x) x$weight))
            mean.weight = mean(pooled.elements.weights)
            pooled.elements.weights = pooled.elements.weights/sum(pooled.elements.weights)
            pooled.mean = colSums(t(pooled.elements.weighted.means)*pooled.elements.weights)
            template.cell = list(mean = pooled.mean, cov = optimalFlow::GaussianBarycenters(pooled.elements.covs,
                                                                               pooled.elements.weights)$Barycenter, weight = mean.weight, type = k, type.proportions =
                                   pooled.elements.type.proportion)
          }
          parallel::stopCluster(cl)
          templates[[j]] = template
          templates.cell.types[[j]] = template.partition.labels
        }
      }
    } else{
      j = 0
      for (i in sort(unique(cytos.cluster))){
        if(i == 0){
          next
        } else{
          j = j + 1
          if(sum(cytos.cluster == i) == 1){
            cell.types = lapply(database.elliptical[cytos.cluster == i], function(x) lapply(x, function(y) y$type))
            cell.types = sort(unique(unlist(cell.types)))
            template.partition = cell.types
            template.partition.labels = sort(unique(template.partition))
            if( length(which(template.partition.labels == 0)) > 0){
              template.partition.labels = template.partition.labels[-1]
            }
            cl <- parallel::makePSOCKcluster(cl.paral)
            doParallel::registerDoParallel(cl)
            template = foreach::foreach(k = template.partition.labels) %dopar%{
              # for(k in template.partition.labels){
              pooled.elements = pooled.clustered.cytometries[[as.character(i)]][template.partition == k]
              pooled.elements.type.proportion = unlist(lapply(pooled.elements, function(x) x$type))
              pooled.elements.type.proportion = table(pooled.elements.type.proportion)
              pooled.elements.type.proportion = sort(pooled.elements.type.proportion/ sum(pooled.elements.type.proportion), decreasing = TRUE)
              pooled.elements.covs = lapply(pooled.elements, function(x) x$cov)
              pooled.elements.weighted.means = matrix(unlist(lapply(pooled.elements, function(x) x$mean)),
                                                      ncol = length(pooled.elements), byrow = FALSE)
              pooled.elements.weights = unlist(lapply(pooled.elements, function(x) x$weight))
              mean.weight = mean(pooled.elements.weights)
              pooled.elements.weights = pooled.elements.weights/sum(pooled.elements.weights)
              pooled.mean = colSums(t(pooled.elements.weighted.means)*pooled.elements.weights)
              template.cell = list(mean = pooled.mean, cov = optimalFlow::GaussianBarycenters(pooled.elements.covs,
                                                                                 pooled.elements.weights)$Barycenter, weight = mean.weight, type = k, type.proportions =
                                     pooled.elements.type.proportion)
            }
            parallel::stopCluster(cl)
            templates[[j]] = template
            templates.cell.types[[j]] = template.partition.labels
          } else{
            if(bar.repetitions < cl.paral){
              n.cores = bar.repetitions
            } else{
              n.cores = cl.paral
            }
            print(length(pooled.clustered.cytometries[[as.character(i)]]))
            if(op.syst == "unix"){
              k.barycenter.list = parallel::mclapply(1:bar.repetitions, optimalFlow::wassersteinKBarycenter, k = barycenters.number, alpha = alpha.bar,
                                                     initialization = "plus-plus", pooled.clustered.cytometries[[as.character(i)]],
                                                     mc.cores = n.cores)#n.cores
            } else{
              cl <- parallel::makePSOCKcluster(n.cores)
              k.barycenter.list = parallel::parLapply(cl, 1:bar.repetitions, optimalFlow::wassersteinKBarycenter, k = barycenters.number, alpha = alpha.bar,
                                                     initialization = "plus-plus", pooled.clustered.cytometries[[as.character(i)]])
              parallel::stopCluster(cl)
            }
            wasser.var = lapply(1:bar.repetitions, function (ii) k.barycenter.list[[ii]]$wasserstein.var)
            if(min(unlist(wasser.var)) == Inf){
              stop("All executions of k-barycenter were unsuccesfull")
            } else{
              print(unlist(wasser.var))
              k.barycenter = k.barycenter.list[[which.min(wasser.var)]]$wasserstein.k.barycenter
              template.partition = k.barycenter$cluster
              template.partition.labels = sort(unique(template.partition))
              print(template.partition)
              if( length(which(template.partition.labels == 0)) > 0){
                template.partition.labels = template.partition.labels[-1]
              }

              if(length(template.partition.labels) < cl.paral){
                n.cores = length(template.partition.labels)
              } else{
                n.cores = cl.paral
              }
              if(op.syst == "unix"){
                template = parallel::mclapply(template.partition.labels, function(k){
                  pooled.elements = pooled.clustered.cytometries[[as.character(i)]][template.partition == k]
                  pooled.elements.type.proportion = unlist(lapply(pooled.elements, function(x) x$type))
                  pooled.elements.type.proportion = table(pooled.elements.type.proportion)
                  pooled.elements.type.proportion = sort(pooled.elements.type.proportion/ sum(pooled.elements.type.proportion), decreasing = TRUE)
                  pooled.elements.weights = unlist(lapply(pooled.elements, function(x) x$weight))
                  mean.weight = mean(pooled.elements.weights)
                  template.cell = list(mean = k.barycenter$baricentro[[which(template.partition.labels == k)]]$mean,
                                       cov = k.barycenter$baricentro[[which(template.partition.labels == k)]]$cov, weight = mean.weight,
                                       type = k, type.proportions = pooled.elements.type.proportion)
                },
                mc.cores = n.cores)
              } else{
                cl <- parallel::makePSOCKcluster(n.cores)
                template = parallel::parLapply(cl, template.partition.labels, function(k){
                  pooled.elements = pooled.clustered.cytometries[[as.character(i)]][template.partition == k]
                  pooled.elements.type.proportion = unlist(lapply(pooled.elements, function(x) x$type))
                  pooled.elements.type.proportion = table(pooled.elements.type.proportion)
                  pooled.elements.type.proportion = sort(pooled.elements.type.proportion/ sum(pooled.elements.type.proportion), decreasing = TRUE)
                  pooled.elements.weights = unlist(lapply(pooled.elements, function(x) x$weight))
                  mean.weight = mean(pooled.elements.weights)
                  template.cell = list(mean = k.barycenter$baricentro[[which(template.partition.labels == k)]]$mean,
                                       cov = k.barycenter$baricentro[[which(template.partition.labels == k)]]$cov, weight = mean.weight,
                                       type = k, type.proportions = pooled.elements.type.proportion)
                })
                parallel::stopCluster(cl)
              }
              templates[[j]] = template
              templates.cell.types[[j]] = template.partition.labels
            }
          }
        }
      }
    }
  }
  sys.time.step2.1 = Sys.time()
  time.dif.2 =  difftime(sys.time.step2.1, sys.time.step2.0)
  print(paste("step 2:", time.dif.2, units(time.dif.2), sep = " "))
  t1 = Sys.time()
  time.dif.total = difftime(sys.time.step2.1, t0)
  print(paste("Execution time:", time.dif.total, units(time.dif.total), sep = " "))
  return(list(templates = templates, clustering = cytos.cluster, database.elliptical = database.elliptical))
}
