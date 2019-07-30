qdaClassification = function(normal, data){
  data = as.matrix(data)
  media = matrix(normal$mean, nrow = 1, ncol = length(normal$mean))
  inv = tryCatch(solve(normal$cov), error = function(x) NA)
  if(is.na(inv)){
    return(inv)
  } else{
    weight = normal$weight
    qda_score = -((1/2)*rowSums((data%*%inv)*data)) + data%*%inv%*%t(media) + matrix(-(1/2)*media%*%inv%*%t(media) -
                (1/2)*log(1/det(inv))+log(weight), nrow = dim(data)[1], ncol = 1)

  }
}
