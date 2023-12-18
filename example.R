
library(Escort)
library(parallelDist) # for fast calculation of the distance matrix.
# library(cluster)
library(mclust)
# library(RMTstat)



# setwd("/Users/ruby/Desktop/Graduate/Research/Escort/")
data("step0_clean_dyntoy_L3")


#########################
######### step1 #########
#########################

# test distinct clusters
dist_mat <- parallelDist::parDist(t(norm_counts), method = "manhattan")
LvsC <- HD_DCClusterscheck(dist_mat=dist_mat, rawcounts=rawcounts)
LvsC$DCcheck

# test Homogeneous cells
cor_test <- step1_testHomogeneous(norm_counts=norm_counts, num.sim = 1000)
cor_test$decision


#########################
######### step2 #########
#########################

gene.var <- quick_model_gene_var(norm_counts)
genes.HVGs <- rownames(gene.var)[1:1000] # top 1000 genes (approx. 20%)
dimred <- getDR_2D(norm_counts[genes.HVGs,], "PCA")
head(dimred)

# check if embedding has distinct clusters
dist_mat <- dist(dimred, method = "euclidean")
DRLvsC <- LD_DCClusterscheck(dist_mat=dist_mat, DRdims=dimred)
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
cls <- mclust::Mclust(dimred)$classification
ti_out <- slingshot(data=dimred, clusterLabels=cls)
rawpse <- slingPseudotime(ti_out, na=T)
ls_fitLine <- lapply(slingCurves(ti_out), function(x) x$s[x$ord,])

library(grDevices)
brewCOLS <- c("#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#E6F598", "#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2")
colors <- colorRampPalette(brewCOLS)(100)
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
calcScore(scoredf)

