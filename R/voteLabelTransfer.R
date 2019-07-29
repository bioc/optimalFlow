voteLabelTransfer = function(type = "points", test.partition, test.cytometry, test.partition.ellipse, training.cytometries,
                             training.cytometries.barycenter,
                             test = 1, op.syst, cl.paral = 1, equal.weights = FALSE){
  if (type == "points"){
    transported.labels = list()
    if(op.syst == "windows"){
      cl = parallel::makePSOCKcluster(cl.paral)
      for (j in test){
        t0 = Sys.time()
        transported.labels[[as.character(j)]] = parallel::parLapply(cl, training.cytometries, optimalFlow::labelTransfer, test.cytometry,
                                                                    test.partition, equal.weights)
        t1 = Sys.time()
        print(paste(t1-t0, units(difftime(t1,t0))))
      }
      parallel::stopCluster(cl)
    } else {
      if (op.syst == "unix"){
        for (j in test){
          t0 = Sys.time()
          transported.labels[[as.character(j)]] = parallel::mclapply(training.cytometries, optimalFlow::labelTransfer, test.cytometry,
                                                                     test.partition, equal.weights, mc.cores = cl.paral)
          t1 = Sys.time()
          print(paste(t1-t0, units(difftime(t1,t0))))
        }

      }
    }
  } else{
    if (type == "ellipses"){
      test.partition = 1:length(test.partition.ellipse)
      training.cytometries = 1
      transported.labels = list()
      for (j in test){
        t0 = Sys.time()
        transported.labels[[as.character(j)]] = lapply(1, optimalFlow::labelTransferEllipse, test.partition.ellipse,
                                                       training.cytometries.barycenter, equal.weights)
        t1 = Sys.time()
        print(paste(t1-t0, units(difftime(t1,t0))))
      }
    }
  }

  vote.transport.labels = list()
  complete.vote.transport.labels = list()
  for(i in as.character(test)){
    transp_lab_cito = list()
    transp_lab_cito_complete = list()
    for(test.cell in unique(transported.labels[[i]][[1]]$to)){
      vote.on.test.cell = c()
      for(j in 1:length(training.cytometries)){
        receiving.cell.indexes = which(transported.labels[[i]][[j]]$to == test.cell)
        sending.cells = transported.labels[[i]][[j]]$from[receiving.cell.indexes]
        sending.cells.vote.proportion = transported.labels[[i]][[j]]$mass[receiving.cell.indexes]/sum(transported.labels[[i]][[j]]$
                                                                                                        mass[receiving.cell.indexes])
        sending.cells.proportions = array(dim = length(receiving.cell.indexes))
        t = 0
        for (jj in receiving.cell.indexes){
          t = t + 1
          receiving.cell.indexes_jj = which(transported.labels[[i]][[j]]$from == transported.labels[[i]][[j]]$from[jj])
          sending.cells.proportions[t] = transported.labels[[i]][[j]]$mass[jj]/sum(transported.labels[[i]][[j]]$
                                                                                     mass[receiving.cell.indexes_jj])
        }
        # sending.cells_orig_prop = sending.cells.proportions
        vote.on.test.cell = rbind(vote.on.test.cell, cbind(rep(j,t), sending.cells, sending.cells.vote.proportion,
                                                           sending.cells.proportions))
      }
      vote.on.test.cell = data.frame("id" = as.integer(unlist(vote.on.test.cell[,1])), "cell" = vote.on.test.cell[,2],
                                     "vote.proportion" = as.double(unlist(vote.on.test.cell[,3])), "original.proportion" =
                                       as.double(unlist(vote.on.test.cell[, 4])))
      vote.on.test.cell$compound.proportion = vote.on.test.cell$vote.proportion*vote.on.test.cell$original.proportion
      vote.on.test.cell.bis = as.data.frame(vote.on.test.cell %>% group_by(cell) %>% summarise(compound.proportion =
                                                                                                 sum(compound.proportion), simple.proportion = sum(vote.proportion)))
      vote.on.test.cell.bis$compound.proportion = vote.on.test.cell.bis$compound.proportion/sum(vote.on.test.cell.bis$
                                                                                                  compound.proportion)
      vote.on.test.cell.bis$simple.proportion = vote.on.test.cell.bis$simple.proportion/sum(vote.on.test.cell.bis$
                                                                                              simple.proportion)
      transp_lab_cito[[test.cell]] = vote.on.test.cell.bis[order(-vote.on.test.cell.bis$compound.proportion),]
      transp_lab_cito_complete[[test.cell]] = vote.on.test.cell
    }
    vote.transport.labels[[i]] = transp_lab_cito
    complete.vote.transport.labels[[i]] = transp_lab_cito_complete
  }

  # for (j in as.character(test)){
  #   matrix.clas = array(0, dim = c(length(cell.types), length(names(vote.transport.labels[[j]]))))
  #   t = 0
  #   for(i in names(vote.transport.labels[[j]])){
  #     t = t + 1
  #     rest.of.sending.cells = which(cell.types %in% vote.transport.labels[[j]][[i]]$cell)
  #     if(length(rest.of.sending.cells) == 0){
  #       rest.of.sending.cells = cell.types
  #     }
  #     for (k in rest.of.sending.cells){
  #       matrix.clas[k, t] = vote.transport.labels[[j]][[i]]$compound.proportion[vote.transport.labels[[j]][[i]]$cell ==
  #                                                                               cell.types[k]]
  #     }
  #     colnames(matrix.clas) = names(vote.transport.labels[[j]])
  #     rownames(matrix.clas) = cell.types
  #   }
  #   gplots::heatmap.2(matrix.clas, col = colorRampPalette(RColorBrewer::brewer.pal(8, "Blues"))(50), key = T, keysize = 1.1,
  #                     density.info =  "none",trace = "none", dendrogram = "none", Rowv = F, Colv = F, cexRow = 0.65,
  #                     lwid = c(0.4, 2), key.title = NA, key.xlab = NA, lhei = c(0.25,2), margins = c(3.2,10.2),symbreaks = F)
  #
  # }
  return(list(final.vote = vote.transport.labels, complete.vote = complete.vote.transport.labels))
}
