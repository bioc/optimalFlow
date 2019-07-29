labelTransfer = function(training.cytometry, test.cytometry, test.partition, equal.weights = FALSE){
  dim.cyto = dim(test.cytometry)[2]
  training.partition = training.cytometry[, dim(training.cytometry)[2]]
  training.cells = names(table(training.partition))
  M = length(training.cells)
  test.cells = names(table(test.partition))
  N = length(test.cells)
  cost.matrix = array(dim = c(M, N))
  t0 = Sys.time()
  j = 0
  for (cell in training.cells) {
    j = j + 1
    n_1 = sum(training.partition == cell)
    jj = 0
    for (k_m in test.cells){
      jj = jj + 1
      n_2 = sum(test.partition == k_m)
      if(n_1 <= 10000){
        if (n_2 <= 10000){
          cost = sum(Rfast::dista(training.cytometry[training.partition == cell, 1:dim.cyto], test.cytometry[test.partition == k_m, 1:dim.cyto]))/
            (n_1*n_2)
        } else{
          cost = sum(Rfast::dista(training.cytometry[training.partition == cell, 1:dim.cyto], test.cytometry[sample(which(test.partition == k_m), 10000), 1:dim.cyto]))/
            (n_1*10000)
        }
      }else{
        if(n_2 <= 10000){
          cost = sum(Rfast::dista(training.cytometry[sample(which(training.partition == cell), 10000), 1:dim.cyto], test.cytometry[test.partition == k_m, 1:dim.cyto]))/
            (n_2*10000)
        } else{
          cost = sum(Rfast::dista(training.cytometry[sample(which(training.partition == cell), 10000), 1:dim.cyto], test.cytometry[sample(which(test.partition == k_m), 10000), 1:dim.cyto]))/
            (10000*10000)
        }
      }
      cost.matrix[j, jj] = cost
    }
    # print(cost.matrix[j,])
  }

  if(equal.weights){
    A = matrix(1/M, nrow = 1, ncol = M)
    B = matrix(1/N, nrow = 1, ncol = N)
    names(A) = training.cells
    names(B) = test.cells

  } else{
    A = table(training.partition)/sum(table(training.partition))
    B = table(test.partition)/sum(table(test.partition))
  }

  optimal.transp.from.A.to.B = transport::transport(a = A, b = B, costm = cost.matrix)
  optimal.transp.from.A.to.B$to = names(B)[optimal.transp.from.A.to.B$to]
  optimal.transp.from.A.to.B$from = names(A)[optimal.transp.from.A.to.B$from]
  return(optimal.transp.from.A.to.B)
}
