estimationCellBarycenter = function(cell, cytometries){
  inexes.cell = lapply(cytometries, function(x) which(unlist(lapply(x, function(y) y$type)) == cell))
  cell.representatives = list()
  j = 0
  for (i in 1:length(inexes.cell)){
    if (length(inexes.cell[[i]]) == 0 ){
      next
    } else{
      j = j + 1
      cell.representatives[[j]] = cytometries[[i]][[inexes.cell[[i]]]]
    }
  }
  # cell.representatives = lapply(cell.representatives, function(i) lapply(cytometries, function(x) x[[i]]))
  if( length(cell.representatives) > 0){
    weight.cell = unlist(lapply(cell.representatives, function(x) x$weight))
    cov.cell = lapply(cell.representatives, function(x) x$cov)
    cov.barycenter.cell = optimalFlow::GaussianBarycenters(cov.cell, weight.cell/sum(weight.cell))$Barycenter
    mean.cell = lapply(cell.representatives, function(x) x$mean)
    mean.barycenter.cell = colSums(matrix(unlist(lapply(1:length(cell.representatives),
                                                        function(i) mean.cell[[i]]*weight.cell[i]/sum(weight.cell))),
                                          nrow = length(cell.representatives), byrow = TRUE))
    weight.barycenter.cell = mean(weight.cell)
    return(list(mean = mean.barycenter.cell, cov = cov.barycenter.cell, weight = weight.barycenter.cell, type = cell))
  }
}
