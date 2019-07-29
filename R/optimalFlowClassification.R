optimalFlowClassification <- function(X, database, templates, consensus.method = "pooling", cov.estimation = "standard", alpha.cov = 0.85,
                                      initial.method = "supervized", alpha.tclust = 0, restr.factor.tclust = 1000,
                                      classif.method = "qda", qda.bar = TRUE, cost.function = "points", cl.paral = 1,
                                      equal.weights.voting = TRUE, equal.weights.template = TRUE){
  t0 = Sys.time()
  op.syst = .Platform$OS.type
  #############################################################################################################
  ############################# Step 0: Checking that entries are well defined ################################
  #############################################################################################################

  if(!is.matrix(X) & !is.data.frame(X)){
    stop("Data is not a matrix or a dataframe.")
  }
  cov.estimation = match.arg(cov.estimation, c("standard", "robust"))
  if(!is.double(alpha.cov)){
    stop("alha.cov is not a doble")
  } else{
    if(alpha.cov < 0 | alpha.cov >= 1){
      stop("alpha.cov is not in the range (0,1).")
    }
  }
  initial.method = match.arg(initial.method, c("supervized", "unsupervized"))
  if(!is.double(alpha.tclust)){
    stop("alha.tclust is not a doble")
  } else{
    if(alpha.tclust < 0 | alpha.tclust >= 1){
      stop("alpha.tclust is not in the range [0,1).")
    }
  }
  classif.method = match.arg(classif.method, c("matching", "qda", "random forest"))
  cost.function = match.arg(cost.function, c("points", "ellipses"))
  op.syst = match.arg(op.syst, c("unix", "windows"))

  #############################################################################################################
  ######################################### Step 1: tclust step ###############################################
  #############################################################################################################
  sys.time.step1.0 = Sys.time()
  dim.cyto = dim(X)[2]
  if(initial.method == "unsupervized"){
    stop("unsupervized initialization is not yet implemented.")
  } else{
    if(length(templates) < cl.paral){
      n.cores = length(templates)
    } else{
      n.cores = cl.paral
    }
    if (op.syst == "unix"){
      tclust.results = parallel::mclapply(templates,
                                function(template, alpha.tclust, restr.factor.tclust){
                                  initial.weights = unlist(lapply(template, function(x) x$weight))
                                  initial.means = matrix(unlist(lapply(template, function(x) x$mean)), ncol = length(initial.weights), byrow = FALSE)

                                  initial.covs = array(dim = c(dim.cyto, dim.cyto, length(initial.weights)))
                                  for( i in 1:length(initial.weights)){
                                    initial.covs[,,i] = template[[i]]$cov
                                  }

                                  template.tclust = list()
                                  template.tclust$cov = initial.covs
                                  template.tclust$centers = initial.means
                                  template.tclust$weights = initial.weights/sum(initial.weights)
                                  clustering.tclust = optimalFlow::tclustWithInitialization(template.tclust, X, i.sol.type = "barycenter",
                                                                               trimming = alpha.tclust, restr.fact = restr.factor.tclust)
                                  return(clustering.tclust)
                                },
                                alpha.tclust = alpha.tclust, restr.factor.tclust = restr.factor.tclust, mc.cores = n.cores)
    } else{
      cl <- parallel::makePSOCKcluster(n.cores)
      tclust.results = parallel::parLapply(cl, templates,
                                function(template, alpha.tclust, restr.factor.tclust){
                                  initial.weights = unlist(lapply(template, function(x) x$weight))
                                  initial.means = matrix(unlist(lapply(template, function(x) x$mean)), ncol = length(initial.weights), byrow = FALSE)

                                  initial.covs = array(dim = c(dim.cyto, dim.cyto, length(initial.weights)))
                                  for( i in 1:length(initial.weights)){
                                    initial.covs[,,i] = template[[i]]$cov
                                  }

                                  template.tclust = list()
                                  template.tclust$cov = initial.covs
                                  template.tclust$centers = initial.means
                                  template.tclust$weights = initial.weights/sum(initial.weights)
                                  clustering.tclust = optimalFlow::tclustWithInitialization(template.tclust, X, i.sol.type = "barycenter",
                                                                               trimming = alpha.tclust, restr.fact = restr.factor.tclust)
                                  return(clustering.tclust)
                                },
                                alpha.tclust = alpha.tclust, restr.factor.tclust = restr.factor.tclust)
      parallel::stopCluster(cl)
    }
    objective.tclust = unlist(lapply(tclust.results, function(x) x$obj))
    assigned.tclust.index = which.max(objective.tclust)

  }

  sys.time.step1.1 = Sys.time()
  time.dif.1 = difftime(sys.time.step1.1, sys.time.step1.0)
  print(paste("step 1:", time.dif.1, units(time.dif.1), sep = " "))
  #############################################################################################################
  ############################### Step 2: Similarity distance assignation #####################################
  #############################################################################################################
  sys.time.step2.0 = Sys.time()

  data.elliptical = lapply(sort(unique(tclust.results[[assigned.tclust.index]]$cluster)), optimalFlow::estimCovCellGeneral, cytometry =
                             X, labels = tclust.results[[assigned.tclust.index]]$cluster, type = cov.estimation,
                           alpha = alpha.cov)
  data.elliptical = data.elliptical[!is.na(data.elliptical)]

  similarity.distances = optimalFlow::costWasserMatchingEllipse(data.elliptical, templates, equal.weights.template)
  print(similarity.distances)
  assigned.template.index = which.min(similarity.distances)
  assigned.template = templates[[assigned.template.index]]

  sys.time.step2.1 = Sys.time()
  time.dif.2 = difftime(sys.time.step2.1, sys.time.step2.0)
  print(paste("step 2:", time.dif.2, units(time.dif.2), sep = " "))
  #############################################################################################################
  ############################################ Step 3: Labelling ##############################################
  #############################################################################################################
  sys.time.step3.0 = Sys.time()

  if(classif.method == "qda" & qda.bar){
    weight.groups = unlist(lapply(assigned.template, function(x) x$weight))
    if(abs(sum(weight.groups) - 1) > 10^(-8)){
      for(i in 1:length(assigned.template)){
        assigned.template[[i]]$weight = assigned.template[[i]]$weight/sum(weight.groups)
      }
    }
    if(length(assigned.template) < cl.paral){
      n.cores = length(assigned.template)
    } else{
      n.cores = cl.paral
    }
    if( op.syst == "unix"){
      qda.assignation = parallel::mclapply(assigned.template, optimalFlow::qdaClassification, data = X,
                                           mc.cores = n.cores)
    } else{
      cl <- parallel::makePSOCKcluster(n.cores)
      qda.assignation = parallel::parLapply(cl, assigned.template, optimalFlow::qdaClassification, data = X)
      parallel::stopCluster(cl)
    }
    qda.assignation = matrix(unlist(qda.assignation), ncol = length(assigned.template), byrow = FALSE)
    qda.groups = Rfast::rowMaxs(qda.assignation)
    unique.qda.groups = sort(unique(qda.groups))
    qda.groups = factor(qda.groups)
    levels(qda.groups) = unlist(lapply(assigned.template, function(x) x$type))[unique.qda.groups]

    if (consensus.method == "hierarchical" |  consensus.method == "k-barycenter"){
      vote = list()
      if(sum(is.na(as.integer(levels(qda.groups)))) > 0){
        range.qda.groups = levels(qda.groups)
      } else{
        range.qda.groups = as.integer(levels(qda.groups))
      }
      for (i in range.qda.groups){
        vote[[as.character(i)]] = data.frame(cell = names(assigned.template[[which(range.qda.groups == i)]]$type.proportions),
                                             simple.proportion =
                                               as.vector(assigned.template[[which(range.qda.groups == i)]]$type.proportions))
      }
      sys.time.step3.1 = Sys.time()
      time.dif.3 = difftime(sys.time.step3.1, sys.time.step3.0)
      print(paste("step 3:", time.dif.3, units(time.dif.3), sep = " "))
      t1 = Sys.time()
      print(difftime(t1, t0))
      return(list(cluster = qda.groups, cluster.vote = vote, vote.original.proportion = lapply(assigned.template, function(x)
        x$type.proportions), clusterings = tclust.results, assigned.template.index = assigned.template.index))
    } else{
      sys.time.step3.1 = Sys.time()
      time.dif.3 = difftime(sys.time.step3.1, sys.time.step3.0)
      print(paste("step 3:", time.dif.3, units(time.dif.3), sep = " "))
      t1 = Sys.time()
      print(difftime(t1, t0))
      return(list(cluster = qda.groups, clusterings = tclust.results, assigned.template.index = assigned.template.index))
    }
  } else{
    if(classif.method == "matching"){
      if(cost.function == "points"){
        cytos.cluster = templates$clustering
        if(sum(cytos.cluster == assigned.template.index) < cl.paral){
          n.cores = sum(cytos.cluster == assigned.template.index)
        } else{
          n.cores = cl.paral
        }
        vote = optimalFlow::voteLabelTransfer(type = cost.function, test.partition = tclust.results[[assigned.tclust.index]]$cluster,
                                 test.cytometry = X, training.cytometries =
                                   database[cytos.cluster == assigned.template.index], test = 1, op.syst = op.syst,
                                 cl.paral = n.cores, equal.weights = equal.weights.voting)
        sys.time.step3.1 = Sys.time()
        time.dif.3 = difftime(sys.time.step3.1, sys.time.step3.0)
        print(paste("step 3:", time.dif.3, units(time.dif.3), sep = " "))
        t1 = Sys.time()
        print(difftime(t1, t0))
        return(list(cluster = tclust.results[[assigned.tclust.index]]$cluster, cluster.vote = vote$final.vote[[1]], clusterings = tclust.results,
                     assigned.template.index = assigned.template.index))

      } else{

        vote = optimalFlow::voteLabelTransfer(type = cost.function, test.partition.ellipse = data.elliptical, training.cytometries.barycenter =
                                   assigned.template, test = 1, op.syst = op.syst,
                                   cl.paral = 1, equal.weights = equal.weights.voting)
        sys.time.step3.1 = Sys.time()
        time.dif.3 = difftime(sys.time.step3.1, sys.time.step3.0)
        print(paste("step 3:", time.dif.3, units(time.dif.3), sep = " "))
        t1 = Sys.time()
        print(difftime(t1, t0))
        if (consensus.method == "hierarchical" |  consensus.method == "k-barycenter"){
          vote.original.proportion = lapply(assigned.template, function(x) x$type.proportions)
          if(sum(cytos.cluster == assigned.template.index) == 1){
            # print(vote)
            # print(vote$final.vote[[1]])
            # final.vote = vote$fianl.vote[[1]]
            # print(final.vote)
            return(list(cluster = tclust.results[[assigned.tclust.index]]$cluster, cluster.vote = vote$final.vote[[1]], clusterings =
                          tclust.results, assigned.template.index = assigned.template.index))
          } else{
            final.vote = voteTransformation(vote$final.vote[[1]], vote.original.proportion)
            return(list(cluster = tclust.results[[assigned.tclust.index]]$cluster, cluster.vote = final.vote, clusterings =
                          tclust.results, assigned.template.index = assigned.template.index))
          }
        } else{
          return(list(cluster = tclust.results[[assigned.tclust.index]]$cluster, cluster.vote = vote$final.vote[[1]], clusterings = tclust.results,
                      assigned.template.index = assigned.template.index))
        }
      }
    } else{
      cytos.cluster = templates$clustering
      database.elliptical = templates$database.elliptical
      similarity.distances = optimalFlow::costWasserMatchingEllipse(data.elliptical, database.elliptical[cytos.cluster == assigned.template.index],
                                                       equal.weights.template)
      assigned.training.cytometry.index = which.min(similarity.distances)
      if(classif.method == "qda"){
        assigned.training.cytometry = database.elliptical[cytos.cluster == assigned.template.index][[assigned.training.cytometry.index]]
        if(length(assigned.training.cytometry) < cl.paral){
          n.cores = length(assigned.training.cytometry)
        } else{
          n.cores = cl.paral
        }
        if( op.syst == "unix"){
          qda.assignation = parallel::mclapply(assigned.training.cytometry, optimalFlow::qdaClassification, data = X, mc.cores = n.cores)
        } else{
          cl <- parallel::makePSOCKcluster(n.cores)
          qda.assignation = parallel::parLapply(cl, assigned.training.cytometry, optimalFlow::qdaClassification, data = X)
          parallel::stopCluster(cl)
        }
        qda.assignation = matrix(unlist(qda.assignation), ncol = length(assigned.training.cytometry), byrow = FALSE)
        qda.groups = Rfast::rowMaxs(qda.assignation)
        unique.qda.groups = sort(unique(qda.groups))
        qda.groups = factor(qda.groups)
        levels(qda.groups) = unlist(lapply(assigned.training.cytometry, function(x) x$type))[unique.qda.groups]

        sys.time.step3.1 = Sys.time()
        time.dif.3 = difftime(sys.time.step3.1, sys.time.step3.0)
        print(paste("step 3:", time.dif.3, units(time.dif.3), sep = " "))
        t1 = Sys.time()
        print(difftime(t1, t0))
        return(list(cluster = qda.groups, clusterings = tclust.results, assigned.template.index = assigned.template.index))
      } else{
        assigned.training.cytometry = database[cytos.cluster == assigned.template.index][[assigned.training.cytometry.index]]
        # ytest = as.factor(citometria$`Population ID (name)`)
        dim.cyto = dim(assigned.training.cytometry)[2]
        y = as.factor(dplyr::pull(assigned.training.cytometry, dim.cyto))
        training.rforest = randomForest::randomForest(assigned.training.cytometry[, 1:(dim.cyto - 1)], y = y, keep.forest = TRUE)
        random.forest.groups = predict(training.rforest, X)
        sys.time.step3.1 = Sys.time()
        time.dif.3 = difftime(sys.time.step3.1, sys.time.step3.0)
        print(paste("step 3:", time.dif.3, units(time.dif.3), sep = " "))
        t1 = Sys.time()
        print(difftime(t1, t0))
        return(list(cluster = random.forest.groups, clusterings = tclust.results, assigned.template.index = assigned.template.index))
      }
    }
  }
}
