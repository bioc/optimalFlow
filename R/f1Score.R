f1Score = function(clustering, cytometry, noise.cells){
  dim.cyto = dim(cytometry)[2]
  prediction_qda_19_con_18 = as.character(clustering)
  prediction_qda_19_con_18[prediction_qda_19_con_18 %in% noise.cells] = rep("noise", sum(prediction_qda_19_con_18 %in% noise.cells))

  clustering = prediction_qda_19_con_18
  lev_qda = unique(clustering)
  p = array(dim = c(1,length(lev_qda)))
  r = p
  t = 0
  for (i in lev_qda){
    t = t + 1
    if (i == "noise"){
      p[t] = sum(clustering == i & dplyr::pull(cytometry, dim.cyto) %in% noise.cells)/sum(clustering == i)
      r[t] = sum(clustering == i & dplyr::pull(cytometry, dim.cyto) %in% noise.cells)/sum(dplyr::pull(cytometry, dim.cyto) %in%
                                                                                                     noise.cells)
    } else{
      p[t] = sum(clustering == i & dplyr::pull(cytometry, dim.cyto) %in% i)/sum(clustering == i)
      r[t] = sum(clustering == i & dplyr::pull(cytometry, dim.cyto) %in% i)/sum(dplyr::pull(cytometry, dim.cyto) %in%
                                                                                   i)
    }
  }
  F1_score = 2*p*r/(p+r)
  F1_score_tot = rbind(F1_score, p, r)
  cels_original = names(table(dplyr::pull(cytometry, dim.cyto)))
  cel_no_encontradas = cels_original[!(cels_original %in% prediction_qda_19_con_18) & !(cels_original %in% noise.cells)]
  F1_score_tot = cbind(F1_score_tot, matrix(0, ncol = length(cel_no_encontradas), nrow = 3))
  colnames(F1_score_tot) = c(lev_qda, cel_no_encontradas)
  rownames(F1_score_tot) = c("F1-score", "Precision", "Recall")
  return(F1_score_tot)
}
