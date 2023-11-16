devtools::install("/Users/ruby/Desktop/Graduate/Research/Escort")


library(Escort)
library(parallelDist) # for fast calculation of the distance matrix.
library(cluster)
library(mclust)
library(RMTstat)
library(scran)


# setwd("/Users/ruby/Desktop/Graduate/Research/Escort/")
load("data/step0_clean_dyntoy_L3.RData")


#########################
######### step1 #########
#########################

# test distinct clusters
dist_mat <- parallelDist::parDist(t(norm_counts), method = "manhattan")
LvsC <- HD_DCClusterscheck(dist_mat=dist_mat, rawcounts=rawcounts)
LvsC$DCcheck

# test Homogeneous cells
cor_test <- testHomogeneous(norm_counts=norm_counts)
cor_test$decision


#########################
######### step2 #########
#########################

gene.var <- modelGeneVar(norm_counts)
genes.HVGs <- getTopHVGs(gene.var, prop=0.2)
dimred <- getDR_2D(norm_counts[genes.HVGs,], "PCA")
head(dimred)

# check if embedding has distinct clusters
DRLvsC <- LD_DCClusterscheck(dist_mat=dist(dimred, method = "euclidean"), DRdims=dimred, connectedCells = 1)
DRLvsC$DCcheck
# check the cell-cell relationship is well-preserved
simi_cells <- Similaritycheck(norm_counts=norm_counts, dimred=dimred, Cluters=LvsC)
simi_cells$GoodRate
# check gof
gof_eval <- GOFeval(dimred)
gof_eval$occupiedRate


#########################
######### step3 #########
#########################

# fit a trajectory
library(slingshot)
cl1 <- Mclust(dimred)$classification
ti_out <- slingshot(data=dimred, clusterLabels=cl1)
rawpse <- slingPseudotime(ti_out, na=T)
ls_fitLine <- lapply(slingCurves(ti_out), function(x) x$s[x$ord,])

library(grDevices)
library(RColorBrewer)
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(rawpse, breaks=100)]

plot(dimred, col = plotcol, pch=16, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
lines(SlingshotDataSet(ti_out), lwd=2, col='black')


fitLine <- do.call(rbind, lapply(ls_fitLine, function(x) {
  df_seg <- cbind(x[-nrow(x),],x[-1,])
  colnames(df_seg) <- c("x0", "y0", "x1", "y1")
  return(df_seg)
}))
obj <- prepTraj(dimred, PT=rawpse, fitLine=fitLine)
ushap_eval <- UshapeDetector(obj)
ushap_eval$Ambpct

##################################
######### scoring system #########
##################################

scoredf <- data.frame(DCcheck=DRLvsC$ifConnected, SimiRetain=simi_cells$GoodRate,
                      GOF=gof_eval$occupiedRate, USHAPE=ushap_eval$Ambpct)
final_df <- score_cal(scoredf)

