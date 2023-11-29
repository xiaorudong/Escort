#' Differentiate Clusters vs Lineages in High-Dimensional Data
#'
#' @param dist_mat A distance object that contains the distance between cells from dist() function
#' @param rawcounts A raw count data matrix: row:genes, column:cells
#' @param K The number of clusters if you have prior information. If left NULL, the number of clusters is estimated by scLCA.
#' @param clust.max The maximum number of clusters used in scLCA. The default number is 10.
#' @param cutoff The density, determined using the Jaccard index cutoff, serves as a criterion for defining the connectivity of cells. The default value is 0.3.
#' @param checkcells The number of cells examined per cluster determines the identification of connected clusters. If left NULL, the default value is calculated as the maximum value of the minimum cluster size divided by 5 and 10.
#' @param connectedCells The minimum number of connected cells between clusters used to identify their connection. If left NULL, the default value is 10.
#' @param checksize The number of neighbors to consider for each check cell to determine their connectivity. If left NULL, the value is set to be equal to "checkcells".
#' @param seed The pseudorandom number used to reproduce the output.
#'
#' @importFrom SingleCellExperiment SingleCellExperiment colData
#' @importFrom scLCA myscLCA
#' @importFrom SC3 sc3
#' @import devtools
#'
#' @return A list about the results of disconnected clusters detection.
#' @export

HD_DCClusterscheck <- function(dist_mat, rawcounts,
                               K=NULL, clust.max=10,
                               cutoff=0.3, checkcells=NULL,
                               connectedCells=NULL, checksize=NULL) {


  if (is.null(K)) {
    # https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02622-0
    myclust.res <- scLCA::myscLCA(rawcounts, clust.max=clust.max)
    c_cl <- myclust.res[[1]]
    names(c_cl) <- colnames(rawcounts)
    K <- max(c_cl)
  } else {

    # https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02622-0
    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7444317/
    # create a SingleCellExperiment object
    sce <- SingleCellExperiment::SingleCellExperiment(
      assays = list(
        counts = as.matrix(rawcounts),
        logcounts = as.matrix(log2(rawcounts + 1))),
      rowData = data.frame(feature_symbol = rownames(rawcounts))
    )


    res <- SC3::sc3(sce, ks = K)
    c_cl_raw <- SingleCellExperiment::colData(res)[,1]
    c_cl <- as.numeric(factor(c_cl_raw, levels = sort(unique(c_cl_raw))))
    names(c_cl) <- rownames(SingleCellExperiment::colData(res))
    K <- max(c_cl)
  }



  if (K==1) {
    DCdecision <- "Only one cluster detected. Cells are likely connected, proceed to the homogeneity check step next."
    return(list(DCcheck=DCdecision, K=K, Clusters=c_cl, ifConnected=NA))
  }
  if (K>1) {
    if (any(table(c_cl)<=10)) {
      totSmall <- sum(table(c_cl)  <= 10)
      return(list(DCcheck=paste0("There are ", K ," clusters detected total. However, ", totSmall,
      " clusters have fewer than < 10 cells. Please examine the data for outliers or try setting K manually."),
                  K=K, Clusters=c_cl, ifConnected=NA))
    }
    if (any(table(c_cl)>10)) {
      Escort::BWClusters_Determination(dist_mat=dist_mat, K=K, c_cl=c_cl,
                               cutoff=cutoff, checkcells=checkcells,
                               connectedCells=connectedCells, checksize=checksize)
    }
  }
}


#' Differentiate Clusters vs Lineages in Embedding
#'
#' @param dist_mat A distance object that contains the distance between cells from dist() function
#' @param DRdims A data frame. a 2D embedding from a DR method
#' @param cutoff The density, determined using the Jaccard index cutoff, serves as a criterion for defining the connectivity of cells. The default value is 0.1.
#' @param max.nc The maximum number of clusters used in NbClust. The default number is 5.
#' @param K The number of clusters if you have prior information. If left NULL, the number of clusters is estimated by NbClust.
#' @param checkcells The number of cells examined per cluster determines the identification of connected clusters. If left NULL, the default value is calculated as the maximum value of the minimum cluster size divided by 5 and 10.
#' @param connectedCells The minimum number of connected cells between clusters used to identify their connection. The default value is 1. If left NULL, the default value is 10.
#' @param checksize The number of neighbors to consider for each check cell to determine their connectivity. If left NULL, the value is set to be equal to "checkcells".
#'
#' @import stats
#'
#' @return A list about the results of disconnected clusters detection.
#' @export


LD_DCClusterscheck <- function(dist_mat, DRdims, cutoff=0.1,
                               max.nc=5, K=NULL, checkcells=NULL,
                               connectedCells=1, checksize=NULL) {

  if(is.null(K)) {
    res.nbclust <- Escort::NbClust(DRdims, distance = "euclidean", min.nc = 2, max.nc = max.nc, method = "complete", index ="all")
    K <- length(unique(res.nbclust$Best.partition))
    c_cl <- res.nbclust$Best.partition
  } else {
    res.hc <- hclust(dist_mat, method = "complete")
    c_cl <- cutree(res.hc, k = K)
  }
  if (any(table(c_cl)<=5)) {
    return(list(DCcheck="There are small clusters deteced. Please do clustering again",
                K=K, Clusters=c_cl))
  } else {
    Escort::BWClusters_Determination(dist_mat=dist_mat, K=K, c_cl=c_cl, cutoff=cutoff,
                             checkcells=checkcells, connectedCells=connectedCells,
                             checksize=checksize)
  }
}





#' Precise Calculation of Cluster vs Lineage Differentiation
#'
#' @param dist_mat A distance object that contains the distance between cells from dist() function
#' @param K The number of clusters identified by preview step.
#' @param c_cl A vector of integers indicating the cluster to which each cell is allocated.
#' @param cutoff The density, determined using the Jaccard index cutoff, serves as a criterion for defining the connectivity of cells. The default value is 0.3.
#' @param checkcells The number of cells examined per cluster determines the identification of connected clusters. If left NULL, the default value is calculated as the maximum value of the minimum cluster size divided by 5 and 10.
#' @param connectedCells The minimum number of connected cells between clusters used to identify their connection. If left NULL, the default value is 10.
#' @param checksize The number of neighbors to consider for each check cell to determine their connectivity. If left NULL, the value is set to be equal to "checkcells".
#'
#' @import S4Vectors
#'
#' @return A list about the results of disconnected clusters detection
#' \itemize{
#'   \item DCcheck - the decision sentence
#'   \item ifConnected - a logical output. TRUE represents no disconnected clusters detected.
#'   \item Jaccardsummary - a table showing any two clusters are connected or not.
#'   \item clusterLocation - a table showing which two clusters are connected.
#'   \item K - The number of clusters identified.
#'   \item Clusters - A vector of integers indicating the cluster to which each cell is allocated.
#'   \item ConnectedCells - a list showing the connected cells between clusters.
#'   \item No_ConnectedCells - the number of the connected cells between clusters.
#' }
#' @export

BWClusters_Determination <- function(dist_mat, K, c_cl, cutoff=0.3, checkcells=NULL,
                                     connectedCells=NULL, checksize=NULL) {

  res <- Escort::BWClusters(dist_mat, c_cl)
  JaccardIndex_ls <- res$JaccardIndex
  JaccardIndex_op_ls <- res$signJaccardIndex
  Dists_W <- res$WithinDist
  Dists_B <- res$BetweenDist

  if(is.null(checkcells)) {
    checkcells <- max(round(min(table(c_cl)/5)), 10)
  }
  topcc <-lapply(JaccardIndex_ls, function(x) lapply(x, function(y) head(sort(y, decreasing = T),checkcells)))

  topcc_tls <-lapply(1:length(topcc), function(x) lapply(1:length(topcc[[x]]), function(y) {
    subls <- topcc[[x]][[y]]
    df <- data.frame(value=names(subls))
    df$ToGroup <- rep(names(topcc[[x]])[y], nrow(df))
    df$InGroup <- rep(names(topcc)[x], nrow(df))
    return(df)
  }))

  topcc_tls <- lapply(topcc_tls, function(x) do.call(rbind, x))
  topcc_df <- as.data.frame(do.call(rbind, topcc_tls))

  combnt <- as.data.frame(t(combn(1:K, 2)))
  ggcomparison_ls <- list()
  connectedc_ls <- list()
  dist_matmat <- as.matrix(dist_mat)
  if(is.null(checksize)) {
    checksize <- checkcells
  }
  for (i in 1:nrow(combnt)) {
    findg <- paste("Group", combnt[i,], sep = "_")
    needrow <- topcc_df$ToGroup %in% findg & topcc_df$InGroup %in% findg
    needc <- topcc_df[needrow,]
    ingroupc <- lapply(1:nrow(needc), function(x) {
      a <- names(head(sort(dist_matmat[needc$value[x], names(c_cl)[c_cl==as.numeric(strsplit(needc$InGroup[x], split = "_")[[1]][2])]]), checksize))
      b <- rep(as.numeric(strsplit(needc$InGroup[x], split = "_")[[1]][2]), length(a))
      names(b) <- a
      return(b)
    })

    betweengroupc <- lapply(1:nrow(needc), function(x) {
      a <- names(head(sort(dist_matmat[needc$value[x], names(c_cl)[c_cl==as.numeric(strsplit(as.character(needc$ToGroup[x]), split = "_")[[1]][2])]]), checksize))
      b <- rep(as.numeric(strsplit(as.character(needc$ToGroup[x]), split = "_")[[1]][2]), length(a))
      names(b) <- a
      return(b)
    })

    connectedc <- list()
    for (each in 1:nrow(needc)) {
      sub_ls <- c(ingroupc[[each]], betweengroupc[[each]])
      c_cl_sub <- sub_ls
      dist_mat_sub <- dist_matmat[names(c_cl_sub), names(c_cl_sub)]
      res_eachc <- Escort::BWClusters(dist_mat_sub, c_cl_sub)
      res_JaccardIndex <- sapply(res_eachc$JaccardIndex, cbind)
      sub <- sapply(res_JaccardIndex, function(x) names(x)[which(x>cutoff)])
      connectedc[[each]] <- as.vector(unlist(sub))
    }
    connectedc_ls[[i]] <- unique(unlist(connectedc))
    needc$ConnectedCellID <- unlist(lapply(connectedc, function(x) paste(x, collapse = ",")))
    ggcomparison_ls[[i]] <- needc
  }
  combnt$Connected <- sapply(connectedc_ls, function(x) {
    if(is.null(connectedCells)) {
      # connectedCells=10
      connectedCells=round(min(table(c_cl)/5))
      if(connectedCells<3) connectedCells==3
      if(connectedCells>10) connectedCells==10
    }
    des_check <- ifelse(length(x)>connectedCells, T, F)
    return(des_check)
  })

  check_mat <- combnt[combnt$Connected,c(1,2)]
  check1 <- length(unique(as.vector(unlist(check_mat))))==K
  check2 <- nrow(check_mat)>=K-1
  check3 <- apply(check_mat, 2, function(x) any(duplicated(x)))
  if(sum(check3)==0) {
    check <- check1 & check2
  }
  if(sum(check3)!= 0) {
    dup_V1 <- check_mat$V1[duplicated(check_mat$V1)]
    dup_V2 <- check_mat$V2[duplicated(check_mat$V2)]
    cg_check <- unique(unlist(check_mat[check_mat$V1 %in% dup_V1 | check_mat$V2 %in% dup_V2, ]))
    dg_check <- unique(unlist(check_mat[!(check_mat$V1 %in% dup_V1 | check_mat$V2 %in% dup_V2), ]))

    check4 <- ifelse(sum(duplicated(c(cg_check, dg_check)))==0, length(cg_check), length(unique(c(cg_check, dg_check))))
    check4 <- check4==K

    check <- check1 & check2 & check4
  }

  DCdecision <- ifelse(check,"Congratulations! Escort did not find null spaces between clusters. Proceed to the homogeneity check step next.",
                       "There appears to be null spaces between clusters. We do not recommend proceeding with trajectory analysis without further investigation. Please see the vignette for recommendations.")

  return(list(DCcheck=DCdecision, ifConnected=check, Jaccardsummary=combnt,
              clusterLocation=check_mat, K=K, Clusters=c_cl,
              ConnectedCells=ggcomparison_ls, No_ConnectedCells=connectedc_ls))
}


#' Preliminary Calculation of Cluster vs Lineage Differentiation
#'
#' @param dist_mat A distance object that contains the distance between cells from dist() function
#' @param c_cl A vector of integers indicating the cluster to which each cell is allocated.
#'
#' @return A list about the results of Cluster vs Lineage Differentiation
#' \itemize{
#'   \item JaccardIndex - The Jaccard Index between cells
#'   \item signJaccardIndex - Edited Jaccard Index. The negative sign is added to show the cell is very close to other clusters.
#'   \item WithinDist - 	A distance list containing cell-to-cell distances within clusters
#'   \item BetweenDist - A distance list containing cell-to-cell distances between clusters
#' }
#' @export

BWClusters <- function(dist_mat, c_cl) {

  dist_mat2 <- as.matrix(dist_mat)

  Dists_W <- list()
  Dists_B <- list()
  JaccardIndex_ls <- list()
  JaccardIndex_op_ls <- list()

  for (i in unique(c_cl)) {
    sub_c <- dist_mat2[which(c_cl==i), which(c_cl==i)]
    Dists_W_cells <- lapply(1:nrow(sub_c), function(x) sub_c[x,-x])
    names(Dists_W_cells) <- rownames(sub_c)
    Dists_W[[paste("Group", i, sep="_")]] <- Dists_W_cells

    notsub_c_ls <- list()
    gp <- c(unique(c_cl))[! unique(c_cl) %in% i]

    for (gp1 in gp) {
      notsub_c <- dist_mat2[which(c_cl==i), which(c_cl==gp1)]
      Dists_B_cells <- lapply(1:nrow(notsub_c), function(x) notsub_c[x,])
      names(Dists_B_cells) <- rownames(notsub_c)
      notsub_c_ls[[paste("Group", gp1, sep = "_")]] <- Dists_B_cells
    }
    Dists_B[[paste("Group", i, sep="_")]] <- notsub_c_ls

    JaccardIndex_gp_ls <- list()
    JaccardIndex_gp_op_ls <- list()
    for (gp1 in names((Dists_B[[paste("Group", i, sep="_")]]))) {
      Dists_B_cells <- Dists_B[[paste("Group", i, sep="_")]][[gp1]]
      JaccardIndex_vec <- c()
      JaccardIndex_op_vec <- c()
      for (c in 1:length(Dists_B_cells)) {
        w_vec <- Dists_W_cells[[c]]
        b_vec <- Dists_B_cells[[c]]

        Jaccard <- Escort::JaccardIndex_fun(w_vec, b_vec, plot=F)
        JaccardIndex_vec <- c(JaccardIndex_vec, Jaccard$JaccardIndex)
        JaccardIndex_op_vec <- c(JaccardIndex_op_vec, Jaccard$signJaccardIndex)
      }
      names(JaccardIndex_vec) <- names(Dists_B_cells)
      names(JaccardIndex_op_vec) <- names(Dists_B_cells)
      JaccardIndex_vec <- ifelse(JaccardIndex_vec<0, 1, JaccardIndex_vec)
      JaccardIndex_gp_ls[[gp1]] <- JaccardIndex_vec
      JaccardIndex_gp_op_ls[[gp1]] <- JaccardIndex_op_vec
    }
    JaccardIndex_ls[[paste("Group", i, sep="_")]] <- JaccardIndex_gp_ls
    JaccardIndex_op_ls[[paste("Group", i, sep="_")]] <- JaccardIndex_gp_op_ls
  }
  return(list(JaccardIndex=JaccardIndex_ls, signJaccardIndex=JaccardIndex_op_ls,
              WithinDist=Dists_W, BetweenDist=Dists_B))
}




###
#' Jaccard Index of cells
#'
#' @param w_vec A distance vector containing cell-to-cell distances within clusters
#' @param b_vec A distance vector containing cell-to-cell distances between clusters
#' @param plot a logical parameter. If left NULL, the defaul setting is FALSE and no plot will be printed.
#'
#' @importFrom sfsmisc integrate.xy
#'
#' @return A list containing the Jaccard Index between cells
#' \itemize{
#'   \item JaccardIndex - The Jaccard Index between cells
#'   \item signJaccardIndex - Edited Jaccard Index. The negative sign is added to show the cell is very close to other clusters.
#' }
#' @export

JaccardIndex_fun <- function(w_vec, b_vec, plot=F) {

  lower <- min(c(w_vec, b_vec)) - 0.2
  upper <- max(c(w_vec, b_vec)) + 0.2

  if(length(w_vec)==1 | length(b_vec)==1) {
    JaccardIndex=0
    JaccardIndex_op=0
  } else {
    # generate kernel densities
    dw <- density(w_vec, from=lower, to=upper)
    db <- density(b_vec, from=lower, to=upper)
    d <- data.frame(x=dw$x, a=dw$y, b=db$y)
    if (plot) {
      plot(dw, main="overlaid density plots")
      lines(db, col="blue")
      legend("topleft", legend=c(deparse(substitute(w_vec)), deparse(substitute(b_vec))),
             col=c("black", "blue"), lty=1:1, cex=0.5)
    }

    # calculate intersection densities
    d$w <- pmin(d$a, d$b)

    # integrate areas under curves
    total <- sfsmisc::integrate.xy(d$x, d$a) + sfsmisc::integrate.xy(d$x, d$b)
    intersection <- sfsmisc::integrate.xy(d$x, d$w)

    # compute overlap coefficient
    JaccardIndex <- intersection/(total-intersection)
    JaccardIndex_op <- ifelse(dw$x[which.max(dw$y)]<=db$x[which.max(db$y)], JaccardIndex, -1*JaccardIndex)
  }
  return(list(JaccardIndex=JaccardIndex, signJaccardIndex=JaccardIndex_op))
}


