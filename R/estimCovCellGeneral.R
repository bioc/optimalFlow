estimCovCellGeneral = function(cell, cytometry, labels, type = "standard", alpha = 0.85){
  indices = which(labels == cell)
  ddd = dim(cytometry)
  size = length(indices)
  weight = size/ddd[1]
  cytometry = cytometry[indices, 1:(ddd[2])]
  if(size >= 2){
    if(type == "standard"){
      return(list(mean = colMeans(cytometry), cov = cov(cytometry), weight = weight, type = cell))
    } else {
      if(type == "robust" & alpha*size > (ddd[2] + 1)){
        rob_est = tryCatch(robustbase::covMcd(cytometry, alpha = alpha), error = function(x) NA)
        if(is.na(rob_est)){
          return(NA)
        } else{
          return(list(mean = rob_est$center, cov = rob_est$cov, weight = weight, type = cell))
        }
      } else{
        return(NA)
      }
    }
  } else{
    return(NA)
  }
}