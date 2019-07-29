distGaussianMean <- function(N1,N2){
  di.mean<-sum((N1$mean-N2$mean)^2)
  distancia<-di.mean
  distancia
}

distGaussianCov <- function(N1,N2){
  A<-N1$cov
  B<-N2$cov
  sqrt.matrix<-function(C){
    e <- eigen(C)
    # for(i in 1:length(e$values)){
    #   if(e$values[i] < 0){
    #     e$values[i] = abs(e$values[i])
    #   }
    # }
    e$values = abs(e$values)
    sqrtA<-e$vectors %*% diag(sqrt(e$values)) %*% t(e$vectors)
    sqrtA
  }
  sqrtA<-sqrt.matrix(A)
  M1<-sqrtA%*%B%*%sqrtA
  M2<-sqrt.matrix(M1)
  M3<-A+B-2*M2
  di.cov<-sum(diag(M3))
  distancia<-di.cov
  distancia
}

distGaussian<-function(N1,N2){
  di.mean<-sum((N1$mean-N2$mean)^2)
  A<-N1$cov
  B<-N2$cov
  sqrt.matrix<-function(C){
    e <- eigen(C)
    # for(i in 1:length(e$values)){
    #   if(e$values[i] < 0){
    #     e$values[i] = abs(e$values[i])
    #   }
    # }
    e$values = abs(e$values)
    sqrtA<-e$vectors %*% diag(sqrt(e$values)) %*% t(e$vectors)
    sqrtA
  }
  sqrtA<-sqrt.matrix(A)
  M1<-sqrtA%*%B%*%sqrtA
  M2<-sqrt.matrix(M1)
  M3<-A+B-2*M2
  di.cov<-sum(diag(M3))
  distancia<-di.mean+di.cov
  distancia
}

GaussianBarycenters<-function(matrices, weight){
  d<-ncol(matrices[[1]])
  k<-length(weight)
  Sn<-diag(d)
  trfix<-0
  for (j in 1:k) trfix<-trfix+sum(weight[j]*sum(diag(matrices[[j]])))
  main.iter<-function(S){
    VS<-sum(diag(S))+trfix
    aux1<-eigen(S)
    # for(i in 1:length(aux1$values)){
    #   if(aux1$values[i] < 0){
    #     aux1$values[i] = abs(aux1$values[i])
    #   }
    # }
    aux1$values = abs(aux1$values)
    S12<-aux1$vectors%*%diag(aux1$values^{1/2})%*%t(aux1$vectors)
    S12inv<-aux1$vectors%*%diag(aux1$values^{-1/2})%*%t(aux1$vectors)
    Saux<-list()
    trmen<-0
    for (j in 1:k) {
      temp<-eigen(S12%*%matrices[[j]]%*%S12)
      # for(i in 1:length(temp$values)){
      #   if(temp$values[i] < 0){
      #     temp$values[i] = abs(temp$values[i])
      #   }
      # }
      temp$values = abs(temp$values)
      Saux[[j]]<-temp$vectors%*%diag(temp$values^{1/2})%*%t(temp$vectors)
      Saux[[j]]<-weight[j]*Saux[[j]]
      trj<-sum(diag(Saux[[j]]))
      trmen<-trmen+trj
    }
    VS<-VS-2*trmen
    S<-Reduce('+',Saux)
    list(iterant=S,V=VS)
  }
  nuevo<-main.iter(Sn)
  Vold<-nuevo$V
  Sn<-nuevo$iterant
  tol<-5
  n.iter<-1
  while(tol>1e-9){
    nuevo<-main.iter(Sn)
    n.iter<-n.iter+1
    Vnew<-nuevo$V
    tol<-Vold-Vnew
    Vold<-Vnew
    Sn<-nuevo$iterant
  }
  list(Barycenter=Sn, Variation=Vnew,Num.iter=n.iter)
}

kcenter<-function(points,kk,center.asigned){
  new.bar<-list()
  Trimmed.Variation.Prov<-0
  for (a in 1:kk){
    cluster.a<-points[which(center.asigned==a)]
    n.a<-length(cluster.a)
    means.a<-NULL
    covs.a<-list()
    for (i in 1:n.a){
      means.a<-cbind(means.a,cluster.a[[i]]$mean)
      covs.a[[i]]<-cluster.a[[i]]$cov
    }
    mean.cluster.a<-apply(means.a,1,mean)
    var.means.cluster.a<-mean(apply(means.a^2,2,sum))-sum(mean.cluster.a^2)
    Alg.Baricenter<-GaussianBarycenters(covs.a,weight=rep(1,n.a)/n.a)
    cov.cluster.a<-Alg.Baricenter$Barycenter
    new.bar[[a]]<-list(mean=mean.cluster.a,cov=cov.cluster.a)
    Trimmed.Variation.Prov<-Trimmed.Variation.Prov+Alg.Baricenter$Variation+var.means.cluster.a
  }
  resultado<-list(kcenters=new.bar,t.variation=Trimmed.Variation.Prov)
  resultado
}

wasserMinDist <- function(points,centres)
{
  uu<-length(points)
  kk<-length(centres)
  aux.f1<-function(u){
    N1<-points[[u]]
    aux.f2<-function(j) distGaussian(N1,centres[[j]])
    distancias<-sapply(1:kk,aux.f2)
    resultado<-c(min(distancias),which.min(distancias))
    resultado
  }
  resultado<-sapply(1:uu,aux.f1)
  resultado
}

trimmedMinDist<-function(points, centres, alpha = 0.1){
  asignacion<-wasserMinDist(points,centres)
  orden.dista<-order(asignacion[1,])
  threshold<-sort(asignacion[1,])[floor(length(points)*(1-alpha))]
  asignacion[2,(asignacion[1,]>threshold)]<-0
  resultado<-asignacion[2,]
  resultado
}

