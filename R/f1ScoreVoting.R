f1ScoreVoting = function(voting, clustering, cytometry, nivel_sup, noise.cells){
  correspondance = array(dim = length(voting))
  clusters.to.eliminate = sort(unique(clustering))[-which(sort(unique(clustering)) %in% names(voting))]
  if(length(clusters.to.eliminate) > 0){
    indexes.to.eliminate = which(clustering %in% clusters.to.eliminate)
    clustering = clustering[-indexes.to.eliminate]
    cytometry = cytometry[-indexes.to.eliminate,]
  }
  dim.cyto = dim(cytometry)[2]
  for( i in 1:length(voting)){
    vote_new = voting[[i]]
    noise = which(vote_new$cell %in% noise.cells)
    if(length(noise) > 0){
      if(length(noise) == length(vote_new$cell)){
        correspondance[i] = "noise"
      } else{
        ruido_s = sum(vote_new$compound.proportion[noise])
        vote_new = vote_new[-noise,]
        if (ruido_s > nivel_sup*vote_new[1,2]){
          correspondance[i] = "noise"
        } else{
          if (vote_new[1,2] < nivel_sup*ruido_s){
            correspondance[i] = paste("undecided", vote_new$cell[1], "noise")
          } else{
            if (dim(vote_new)[1]>2){
              if(vote_new[1,2]< nivel_sup*vote_new[2,2]){
                correspondance[i] = paste("undecided", vote_new$cell[1], vote_new$cell[2])
              } else{
                correspondance[i] = paste(vote_new$cell[1])
              }
            } else{
              correspondance[i] = paste(vote_new$cell[1])
            }
          }
        }
      }

    } else{
      if (dim(vote_new)[1]>2){
        if(vote_new[1,2]< nivel_sup*vote_new[2,2]){
          correspondance[i] = paste("undecided", vote_new$cell[1], vote_new$cell[2])
        } else{
          correspondance[i] = paste(vote_new$cell[1])
        }
      } else{
        correspondance[i] = paste(vote_new$cell[1])
      }
    }
  }
  t = 0
  aciertos = array(dim = length(voting))
  for (j in names(voting)){
    t = t + 1
    if (correspondance[t] == "noise"){
      aciertos[t] = sum(dplyr::pull(cytometry, dim.cyto)[clustering == j] %in%
                          noise.cells)/length(dplyr::pull(cytometry, dim.cyto)[clustering == j])
    } else{
      aciertos[t] = sum(dplyr::pull(cytometry, dim.cyto)[clustering == j] %in%
                          correspondance[t])/length(dplyr::pull(cytometry, dim.cyto)[clustering == j])
    }
  }

  p = array(dim = c(1,length(unique(correspondance))))
  r = p
  t = 0
  for (i in unique(correspondance)){
    t = t + 1
    if (i == "noise"){
      p[t] = sum((aciertos*table(clustering))[correspondance %in% "noise"])/sum(table(clustering)[correspondance %in% "noise"])
      r[t] = sum((aciertos*table(clustering))[correspondance %in% "noise"])/sum(dplyr::pull(cytometry, dim.cyto) %in%
                                                                                  noise.cells)
    } else{
      p[t] = sum((aciertos*table(clustering))[correspondance %in% i])/sum(table(clustering)[correspondance %in% i])
      r[t] = sum((aciertos*table(clustering))[correspondance %in% i])/sum(dplyr::pull(cytometry, dim.cyto) %in% i)
    }
  }
  F1_score = 2*p*r/(p+r)
  F1_score_wass = rbind(F1_score, p, r)
  cels_original = names(table(dplyr::pull(cytometry, dim.cyto)))
  cel_no_encontradas = cels_original[!(cels_original %in% correspondance) & !(cels_original %in% noise.cells)]
  F1_score_wass = cbind(F1_score_wass, matrix(0,ncol = length(cel_no_encontradas), nrow = 3))
  colnames(F1_score_wass) = c(unique(correspondance), cel_no_encontradas)
  rownames(F1_score_wass) = c("F1-score", "Precision", "Recall")
  return(list(F1_score = F1_score_wass, correspondance = correspondance))
}
