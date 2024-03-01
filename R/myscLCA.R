#' LCA scRNA analysis
#'
#' This function loads a gene by cell transcript count matarix (row by column). This is from https://bitbucket.org/scLCA/single_cell_lca/src/master/
#'
#' @param datmatrix data matrix, gene by cell
#' @param clust.max maximum number of clusters specified by user, default: 10
#' @param trainingSetSize the number of cells used in training set; default: 1000; maximum: the total number of cells in data matrix
#' @param datBatch a vector of batch group ID for each cell, default: NULL
#' @param zerocorrection correction for zeroes before log transformation, default: 0.25
#' @return an integer vector showing clustering membership
#' @examples
#' # Load example dataset
#' data(myscExampleData)
#' 
#' # Both data matrix and true labels of cells are provided
#' names(myscExampleData)
#' 
#' # With 14,074 genes and 250 cells, 83% of entries are zero
#' dim(myscExampleData$datamatrix)
#' 
#' # Three types of cells
#' table(myscExampleData$truelabel)
#' 
#' # Start my scRNAseq LCA analysis
#' myclust.res <- myscLCA(myscExampleData$datamatrix)
#' 
#' # Top results are provide in a list
#' length(myclust.res)
#' 
#' # The result ranked as the best, compared with the true labels:
#' table(myclust.res[[1]],myscExampleData$truelabel)
#' 
#' 
#' @export 

myscLCA <- function(datmatrix, normalzeddat=NULL, cor.thresh = 0.5, clust.max=10, trainingSetSize=1000, datBatch=NULL, outlier.filter=F, zerocorrection = 0.25){
  myRes <- list()
  
  trainingSetSize = min(ncol(datmatrix),trainingSetSize)
  myDelta = zerocorrection; 
  TPM = datmatrix
  mytmp.sd <- apply(TPM,1,sd)#cc
  TPM = TPM[mytmp.sd>0,]#cc
  UMI = colSums(TPM)
  if(is.null(normalzeddat)){#cc
    TPM.log.all = log2(TPM / matrix(UMI, ncol=ncol(TPM), nrow=nrow(TPM), byrow=T) * 1e4 + myDelta)
  }else{
    TPM.log.all = normalzeddat
  }
  
  UMI.all = UMI
  num.sample = trainingSetSize
  #R.seed=123454321
  #set.seed(R.seed)
  cell.selected = sample(length(UMI.all))[1:num.sample]
  UMI = UMI.all[cell.selected]
  TPM.log = TPM.log.all[, cell.selected]
  mytmp.sd <- apply(TPM.log,1,sd)#cc
  gene.selected = which(mytmp.sd > 0)#cc
  TPM.log = TPM.log[gene.selected,]
  pca = prcomp(TPM.log, scale=T)
  pcaT = prcomp(t(TPM.log), scale=T)
  nTWT = which(p.adjust(TWtest(pcaT$sdev)$pvalue, "BH") <= 0.05)
  if (length(nTWT) < 2) nTWT = 1:2
  nTWT2 = nTWT
  nTWT = nTWT[which(abs(1-cosDist.part(rbind(log2(UMI) - mean(log2(UMI)), t(pcaT$x[,nTWT])), idx1 = 1, idx2 = -1)) < 0.9)]
  if (length(nTWT) == 0) nTWT = max(nTWT2) + 1
  nTW = which(p.adjust(TWtest(pca$sdev)$pvalue, "BH") <= 0.05)
  max.factor = max(nTW)
  factors.best = nTW
  nTW = nTW[which(abs(1-cosDist.part(rbind(log2(UMI), t(pca$rotation[,nTW])), idx1 = 1, idx2 = -1)) < sqrt(cor.thresh))]
  X.rot= pca$rotation
  rots = NULL
  if (is.null(datBatch)) {
    batch1 = matrix(log2(UMI), nrow=1, ncol=ncol(TPM.log))
    #nTW = nTW[which(cor(log2(UMI), pca$rotation[,nTW])^2 < 0.64)]
    ######
  } else {
    batch = datBatch[cell.selected]#sim.data$mycells.ls$batchlabel
    batch1 = matrix(1, nrow= length(unique(batch))+1, ncol = length(batch))
    for (i in 1:length(unique(batch))) batch1[i, which(batch != unique(batch)[i])] = 2
    batch1[nrow(batch1),] = log2(UMI)
  }
  r2.sum = array(0, nrow(batch1))
  for (i in 1:nrow(batch1)) r2.sum[i] = sum(cor(batch1[i,], pca$rotation[, nTW])^2)
  
  batch.order = sort(r2.sum, decreasing=T, index.return=T)$ix
  for (i in batch.order) {
    if (length(nTW) == 0) break
    cor.order = sort(abs(cor(batch1[i,], X.rot[, nTW])), decreasing=T, index.return=T)
    if (length(nTW) == 1) {
      if (cor.order$x[1]^2 < cor.thresh) next
      nTW = NULL; break
    } else {
      if (sum(cor.order$x[1:2]^2) < cor.thresh) next
      for (j in 2:length(nTW)) {
        val2 = cor(batch1[i,], X.rot[, nTW[cor.order$ix[c(1, j)]]])
        val2 = val2/sqrt(sum(val2^2))
        rot=cbind(c(val2[1], val2[2]), c(-val2[2], val2[1]))
        rots = rbind(rots, c(nTW[cor.order$ix[c(1,j)]], val2))
        X.rot[, nTW[cor.order$ix[c(1,j)]]] = X.rot[, nTW[cor.order$ix[c(1,j)]]] %*% rot
      } 
      nTW = nTW[-cor.order$ix[1]]
    }
  }
  
  XT.rot= pcaT$x
  r2.sum = array(0, nrow(batch1))
  for (i in 1:nrow(batch1)) r2.sum[i] = sum(cor(batch1[i,], pcaT$x[, nTWT])^2)
  
  batch.order = sort(r2.sum, decreasing=T, index.return=T)$ix
  for (i in batch.order) {
    if (length(nTWT) == 0) break
    cor.order = sort(abs(cor(batch1[i,], XT.rot[, nTWT])), decreasing=T, index.return=T)
    if (length(nTWT) == 1) {
      if (cor.order$x[1]^2 < cor.thresh) next
      nTWT = NULL; break
    } else {
      if (sum(cor.order$x[1:2]^2) < cor.thresh) next
      for (j in 2:length(nTWT)) {
        val2 = cor(batch1[i,], XT.rot[, nTWT[cor.order$ix[c(1, j)]]])
        val2 = val2/sqrt(sum(val2^2))
        rot=cbind(c(val2[1], val2[2]), c(-val2[2], val2[1]))
        XT.rot[, nTWT[cor.order$ix[c(1,j)]]] = XT.rot[, nTWT[cor.order$ix[c(1,j)]]] %*% rot
      } 
      nTWT = nTWT[-cor.order$ix[1]]
    }
  }
  
  ######
  ######
  #factors.best = nTW
  if (length(nTW) == 0) factors.best = max(factors.best) + 1 else factors.best = nTW
  if (length(nTWT) == 0) nTWT = max(nTWT2) + 1
  if (max.factor < max(factors.best)) max.factor = max(factors.best)
  if (length(nTWT) > 2) {
    pcaT.dis = cosDist(XT.rot[,nTWT]) 
    if (quantile(sqrt(apply(XT.rot[, nTWT]^2, 1, sum)), 0.9) > 4 * median(sqrt(apply(XT.rot[, nTWT]^2, 1, sum)))) pcaT.dise = as.matrix(dist(XT.rot[,nTWT])) else pcaT.dise = NULL
  } else {
    pcaT.dis = as.matrix(dist(XT.rot[,nTWT]))
    pcaT.dise = NULL
  }
  if (length(nTW) > 1) {
    cutoff = quantile(apply(X.rot[,factors.best], 1, sumsq), 1)
    outliers = (apply(X.rot[,factors.best], 1, sumsq) > cutoff)
    if (outlier.filter) rot2 = X.rot[!outliers,factors.best] else rot2 = X.rot[,factors.best]
    myDist.best = cosDist(X.rot[,factors.best])  
    myDist.euc = as.matrix(dist(X.rot[,factors.best]))  
    A = 1 - myDist.best / 2
    D = diag(apply(A, 1, sum))
    L <- (D %^% (-1/2)) %*% A %*% (D %^% (-1/2))
    evL <- eigen(L, symmetric=TRUE)
    nF = clust.max
    if (nF < 10) nF = 10
    if (nF > 50) nF = 50
    if (outlier.filter) X = evL$vectors[!outliers,1:nF] else X = evL$vectors[,1:nF]
    X = X/sqrt(apply(X, 1, sumsq))
    myHClust2 = hclust(dist(X), method="average")
    if (outlier.filter) {
      myHClust = hclust(as.dist(myDist.euc[!outliers, !outliers]), method="average")
    } else {
      myHClust = hclust(as.dist(myDist.euc), method="average")
    }
    sil.best = NA
    sil.best2 = rep(NA, clust.max)
    clusters = matrix(NA, nrow=clust.max, ncol=ncol(TPM.log))  
    for (nClust in 2: clust.max) {
      clust1 = cutree(myHClust, nClust) # inital seed
      mycenter = matrix(NA, nrow=nClust, ncol=length(nTW))
      for (i in 1:nClust) {
        if (sum(clust1 == i) > 1) {
          mycenter[i,] = colMeans(X.rot[which(clust1==i),nTW])
        } else {
          mycenter[i,] = X.rot[which(clust1==i),nTW]
        }
      }
      clust2 = kmeans(X.rot[,nTW], mycenter)$cluster
      clust1 = cutree(myHClust2, nClust)# inital seed
      mycenter = matrix(NA, nrow=nClust, ncol=ncol(X))
      for (i in 1:nClust) {
        if (sum(clust1 == i) > 1) {
          mycenter[i,] = colMeans(X[which(clust1==i),])
        } else {
          mycenter[i,] = X[which(clust1==i),]
        }
      }
      clust3 = kmeans(X, mycenter)$cluster
      if (outlier.filter) {
        sil.1e = mean(cluster::silhouette(clust2, pcaT.dis[!outliers,!outliers])[,3])
        sil.1c = mean(cluster::silhouette(clust3, pcaT.dis[!outliers,!outliers])[,3])
        if (!is.null(pcaT.dise)) {
          sil.c = mean(cluster::silhouette(clust2, pcaT.dise[!outliers,!outliers])[,3])
          if (sil.c > sil.1e) sil.1e = sil.c
          sil.c = mean(cluster::silhouette(clust3, pcaT.dise[!outliers,!outliers])[,3])
          if (sil.c > sil.1c) sil.1c = sil.c
        }
        euc.proj = T
        sil.f = sil.1e
        if (sil.1c > sil.1e) {
          clust2 = clust3
          sil.f = sil.1c
          euc.proj = F
        }
        clust1 = array(NA, nrow(X.rot))
        clust1[!outliers] = clust2
        clust2 = clust1
        for (i in which(outliers)) {
          tmp = array(NA, nClust)
          for (j in 1:nClust)
            if (euc.proj) {
              tmp[j] = mean(myDist.euc[i, which(clust2 == j)])
            } else {
              tmp[j] = mean(myDist.best[i, which(clust2 == j)])
            }
          clust1[i] = which.min(tmp)
        }
      } else {
        sil.1e = mean(cluster::silhouette(clust2, pcaT.dis)[,3])
        sil.1c = mean(cluster::silhouette(clust3, pcaT.dis)[,3])
        if (!is.null(pcaT.dise)) {
          sil.c = mean(cluster::silhouette(clust2, pcaT.dise)[,3])
          if (sil.c > sil.1e) sil.1e = sil.c
          sil.c = mean(cluster::silhouette(clust3, pcaT.dise)[,3])
          if (sil.c > sil.1c) sil.1c = sil.c
        }
        euc.proj = T
        sil.f = sil.1e
        if (sil.1c > sil.1e) {
          clust2 = clust3
          sil.f = sil.1c
          euc.proj = F
        }
        clust1 = clust2
      }
      if (length(unique(clust1)) > 1) {
        sil.best2[nClust] = sil.f
      } else {
        sil.best2[nClust] = 0
      }
      clusters[nClust,] = clust1
      if (is.na(sil.best) | sil.best2[nClust] > sil.best) {
        sil.best = sil.best2[nClust]
        clust.best = clust1
      }
    }
  } else {
    sil.best = NA
    sil.best2 = rep(NA, clust.max)
    clusters = matrix(NA, nrow=clust.max, ncol=ncol(TPM.log))  
    myDist.best = as.matrix(dist(X.rot[,factors.best]))
    cutoff = quantile(X.rot[,factors.best]^2, 1)
    outliers = (X.rot[,factors.best]^2 > cutoff)
    if (outlier.filter) X = X.rot[!outliers,factors.best] else X = X.rot[,factors.best]
    for (nClust in 2:clust.max) {
      clust1 = Mclust(X, G=nClust, verbose = F)$classification
      if (outlier.filter) {
        clust3 = array(NA, ncol(X.rot))
        clust3[!outliers] = clust1
        clust2 = clust3
        for (i in which(outliers)) {
          tmp = array(NA, nClust)
          for (j in 1:nClust) tmp[j] = mean(myDist.best[i, which(clust2 == j)])
          clust3[i] = which.min(tmp)
        }
      } else {
        clust3 = clust1
      }
      if (length(unique(clust3)) > 1) {
        sil.f = mean(cluster::silhouette(clust3, pcaT.dis)[,3])
        if (!is.null(pcaT.dise)) {
          sil.t = mean(cluster::silhouette(clust3, pcaT.dise)[,3])
          if (sil.t > sil.f) sil.f = sil.t
        }
        sil.best2[nClust] = sil.f
      } else {
        sil.best2[nClust] = 0
      }
      clusters[nClust,] = clust3
      if (is.na(sil.best) | sil.best2[nClust] > sil.best) {
        sil.best = sil.best2[nClust]
        clust.best = clust3
      }
    }
  }
  if (length(unique(clust.best)) == 1) {
    myRes[[1]] = rep(1, ncol(TPM))
    return(myRes)
  }
  
  pval = matrix(NA, nrow=length(unique(clust.best)), ncol=max(factors.best))
  for (j in 1:nrow(pval)) for (i in factors.best) pval[j,i] = wilcox.test(X.rot[which(clust.best==unique(clust.best)[j]), i], X.rot[which(clust.best!=unique(clust.best)[j]),  i])$p.value
  factors.best = which(apply(pval, 2, min) < 0.05/(ncol(pval)-1)/nrow(pval))
  nClust = which.max(sil.best2)#nClust = nrow(pval)
  if (length(factors.best) >= 2) {
    cutoff = quantile(apply(X.rot[,factors.best], 1, sumsq), 1)
    outliers = (apply(X.rot[,factors.best], 1, sumsq) > cutoff)
    if (outlier.filter) rot2 = X.rot[!outliers,factors.best] else rot2 = X.rot[,factors.best]
    myDist.best = cosDist(X.rot[,factors.best])
    myDist.euc = as.matrix( dist(X.rot[,factors.best]) ) 
    A = 1 - myDist.best / 2
    D = diag(apply(A, 1, sum))
    L <- (D %^% (-1/2)) %*% A %*% (D %^% (-1/2))
    evL <- eigen(L, symmetric=TRUE)
    nF = length(factors.best)
    if (nF > 50) nF = 50
    if (outlier.filter) X = evL$vectors[!outliers,1:nF] else X = evL$vectors[,1:nF]
    X = X/sqrt(apply(X, 1, sumsq))
    myHClust2 = hclust(dist(X), method="average")
    if (outlier.filter) {
      myHClust = hclust(as.dist(myDist.euc[!outliers, !outliers]), method="average")
    } else {
      myHClust = hclust(as.dist(myDist.euc), method="average")
    }
    clust1 = cutree(myHClust, nClust) # inital seed
    mycenter = matrix(NA, nrow=nClust, ncol=length(nTW))
    for (i in 1:nClust) {
      if (sum(clust1 == i) > 1) {
        mycenter[i,] = colMeans(X.rot[which(clust1==i),nTW])
      } else {
        mycenter[i,] = X.rot[which(clust1==i),nTW]
      }
    }
    clust2 = kmeans(X.rot[,nTW], mycenter)$cluster
    clust1 = cutree(myHClust2, nClust)# inital seed
    mycenter = matrix(NA, nrow=nClust, ncol=ncol(X))
    for (i in 1:nClust) {
      if (sum(clust1 == i) > 1) {
        mycenter[i,] = colMeans(X[which(clust1==i),])
      } else {
        mycenter[i,] = X[which(clust1==i),]
      }
    }
    clust3 = kmeans(X, mycenter)$cluster
    if (outlier.filter) {
      sil.1e = mean(cluster::silhouette(clust2, pcaT.dis[!outliers,!outliers])[,3])
      sil.1c = mean(cluster::silhouette(clust3, pcaT.dis[!outliers,!outliers])[,3])
      if (!is.null(pcaT.dise)) {
        sil.c = mean(cluster::silhouette(clust2, pcaT.dise[!outliers,!outliers])[,3])
        if (sil.c > sil.1e) sil.1e = sil.c
        sil.c = mean(cluster::silhouette(clust3, pcaT.dise[!outliers,!outliers])[,3])
        if (sil.c > sil.1c) sil.1c = sil.c
      }
      euc.proj = T
      sil.f = sil.1e
      if (sil.1c > sil.1e) {
        clust2 = clust3
        sil.f = sil.1c
        euc.proj = F
      }
      clust1 = array(NA, nrow(X.rot))
      clust1[!outliers] = clust2
      clust2 = clust1
      for (i in which(outliers)) {
        tmp = array(NA, nClust)
        for (j in 1:nClust)
          if (euc.proj) {
            tmp[j] = mean(myDist.euc[i, which(clust2 == j)])
          } else {
            tmp[j] = mean(myDist.best[i, which(clust2 == j)])
          }
        clust1[i] = which.min(tmp)
      }
    } else {
      sil.1e = mean(cluster::silhouette(clust2, pcaT.dis)[,3])
      sil.1c = mean(cluster::silhouette(clust3, pcaT.dis)[,3])
      if (!is.null(pcaT.dise)) {
        sil.c = mean(cluster::silhouette(clust2, pcaT.dise)[,3])
        if (sil.c > sil.1e) sil.1e = sil.c
        sil.c = mean(cluster::silhouette(clust3, pcaT.dise)[,3])
        if (sil.c > sil.1c) sil.1c = sil.c
      }
      euc.proj = T
      sil.f = sil.1e
      if (sil.1c > sil.1e) {
        clust2 = clust3
        sil.f = sil.1c
        euc.proj = F
      }
      clust1 = clust2
    }
    clust.best = clust1
  } else if (length(factors.best) == 1) {
    myDist.best = as.matrix(dist(X.rot[,factors.best]))
    if (outlier.filter) X = X.rot[!outliers,factors.best] else X = X.rot[,factors.best]
    clust1 = Mclust(X, G=nClust, verbose = F)$classification
    if (outlier.filter) {
      clust.best = array(NA, ncol(X.rot))
      clust.best[!outliers] = clust1
      clust2 = clust.best
      for (i in which(outliers)) {
        tmp = array(NA, nClust)
        for (j in 1:nClust) tmp[j] = mean(myDist.best[i, which(clust2 == j)])
        clust.best[i] = which.min(tmp)
      }
    } else {
      clust.best = clust1
    }
  } else {
    myRes[[1]] = rep(1, ncol(TPM))
    return(myRes)
  }
  
  if (length(unique(clust.best)) == 1) {
    myRes[[1]] = rep(1, ncol(TPM))
    return(myRes)
  }
  #Projection phase
  mysvd.d=(length(gene.selected)-1)*pca$sdev[1:max.factor]^2
  mytrans.mat = (1/mysvd.d)*t(pca$x[,1:max.factor])
  myCenter = matrix(0, nrow=nClust, ncol = length(factors.best))
  #print(paste("nClust =", nClust, "dim(myCenter) =", dim(myCenter)))
  if (length(factors.best) > 1) {
    for (i in 1:nClust) myCenter[i,] = colMeans(X.rot[which(clust.best == i), factors.best])
  } else {
    for (i in 1:nClust) myCenter[i,] = mean(X.rot[which(clust.best == i), factors.best])
  }
  TPM.log2 = scale(TPM.log.all[gene.selected,])
  myQsize = ncol(TPM.log2)
  myQueries <- t(mytrans.mat %*% TPM.log2)
  if (!is.null(rots)) {
    for (i in 1:nrow(rots)) {
      rot = cbind(c(rots[i,3], rots[i,4]), c(-rots[i,4], rots[i,3]))
      myQueries[, rots[i,1:2]] = myQueries[, rots[i,1:2]] %*% rot
    }
  }
  if (length(factors.best) > 1) {
    if (euc.proj) {
      myassignment = as.matrix(dist(rbind(myQueries[,factors.best],myCenter)))[myQsize+(1:nClust), 1:myQsize]
    } else {
      myassignment = cosDist.part(rbind(myQueries[,factors.best],myCenter),idx1=myQsize+(1:nClust),idx2=1:myQsize)#
    }
  } else {
    myassignment = matrix(nrow=nClust, ncol = myQsize)
    for (i in 1:nClust) myassignment[i,] = abs(myQueries[,factors.best] - as.array(myCenter)[i])
  }
  clust.pred = apply(myassignment, 2, which.min)
  myRes[[1]] <- clust.pred
  return(myRes);
}

