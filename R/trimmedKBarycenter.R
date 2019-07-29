trimmedKBarycenter<-function(k, alpha0, type.ini = "rnd", reps.list){
  ### Primer paso
  if(type.ini == "rnd"){
    indices_prov = sample(1:(length(reps.list)),k, replace = FALSE)
    indices_prov = sort(indices_prov)
    bar.prov<-reps.list[indices_prov]
  }else{
    if(type.ini == "plus-plus"){
      length_rep_lista = length(reps.list)
      j = 1
      c0 = sample(1:(length_rep_lista), 1)
      centers = list(reps.list[[c0]])

      if(k > 1){
        cent_pos = (1:(length_rep_lista))[-c0]
        dist_w2 = array(0, dim = c(k,length_rep_lista))
        dist_w2[1,] = abs(unlist(lapply(reps.list, distGaussian, centers[[1]])))
        D2 = dist_w2[1,][-c0]/sum(dist_w2[1,][-c0])
        # clust = array(dim = u)

        for(j in 2:k){
          # print(D2)
          # print(length(cent_pos))
          c0 = sample(cent_pos, 1, prob = D2)
          c0_pos = which(cent_pos == c0)
          cent_pos = cent_pos[-c0_pos]
          centers[[j]] = reps.list[[c0]]
          dist_w2[j,] = unlist(lapply(reps.list, distGaussian, centers[[j]]))
          s=0
          D2 = array(0, dim = length(cent_pos))
          for(i in cent_pos){
            s=s+1
            # d_min = abs(min(dist_w2[1:j, i]))
            d_min = min(abs(dist_w2[1:j, i]))
            # d_min_clust = which(dist_w2[1:j,i]=d_min)[1]
            D2[s] = d_min
          }
          D2 = D2/sum(D2)
        }
      }
      bar.prov <- centers
    }
    else{
      return("Inicialization not recognised")
    }
    # else{
    #   if (type.ini == "exaus"){
    #     length_rep_lista = length(reps.list)
    #     c0 = l
    #     centers =list(reps.list[[c0]])
    #     if(k > 1){
    #       cent_pos = (1:(length_rep_lista))[-c0]
    #       dist_w2 = array(0, dim = c(k,length_rep_lista))
    #       dist_w2[1,] = unlist(lapply(reps.list, distGaussian, centers[[1]]))
    #       D2 = dist_w2[1,][-c0]/sum(dist_w2[1,][-c0])
    #       # clust = array(dim = u)
    #       for(j in 2:k){
    #         # print(D2)
    #         # print(length(cent_pos))
    #         c0 = sample(cent_pos, 1, prob = D2)
    #         c0_pos = which(cent_pos == c0)
    #         cent_pos = cent_pos[-c0_pos]
    #         centers[[j]] = reps.list[[c0]]
    #         dist_w2[j,] = unlist(lapply(reps.list, distGaussian, centers[[j]]))
    #         s=0
    #         D2 = array(0, dim = length(cent_pos))
    #         for(i in cent_pos){
    #           s=s+1
    #           d_min = abs(min(dist_w2[1:j, i]))
    #           # d_min_clust = which(dist_w2[1:j,i]=d_min)[1]
    #           D2[s] = d_min
    #         }
    #         D2 = D2/sum(D2)
    #       }
    #     }
    #     bar.prov <- centers
    #   }
    # }
  }


  c.asig<-trimmedMinDist(reps.list, bar.prov,alpha= alpha0)
  update.step<-kcenter(reps.list,k,c.asig)
  new.bar<-update.step$kcenters
  Trimmed.Variation.Prov<-update.step$t.variation
  Trimmed.Variation.Old<-Trimmed.Variation.Prov
  bar.prov = new.bar

  ### Iteracion hasta convergencia

  diferencia<-10
  iter<-0
  while(diferencia>0.0000001){
    iter<-iter+1
    c.asig<-trimmedMinDist(reps.list,bar.prov,alpha = alpha0)
    update.step<-kcenter(reps.list,k,c.asig)
    new.bar<-update.step$kcenters
    Trimmed.Variation.Prov<-update.step$t.variation
    diferencia<-Trimmed.Variation.Old-Trimmed.Variation.Prov
    Trimmed.Variation.Old<-Trimmed.Variation.Prov
    bar.prov = new.bar
  }
  Trimmed.Variation<-Trimmed.Variation.Prov
  list(variacion_wasser = Trimmed.Variation, baricentro = bar.prov, cluster = c.asig)
}
