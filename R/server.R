library(parallelDist)
library(mclust)
library(scales)
library(slingshot)
library(DT)
library(RColorBrewer)
library(shinyjs)
library(dplyr)
library(Seurat)
library(clusterProfiler)


server <- function(input, output) {
  # add image:
  output$home_img <- renderImage({

    list(src = "vignettes/Workflow.png",
         width = 600,
         height = 500)
  }, deleteFile = F)

  # step0:
  rawcounts <- reactive({
    req(input$rawfile)
    mat <- read.csv(input$rawfile$datapath, row.names=1)
    as.matrix(mat)
  })

  norm_counts <- reactive({
    req(input$normfile)
    mat <- read.csv(input$normfile$datapath, row.names=1)
    as.matrix(mat)
  })

  # begin analysis
  mydf <- reactiveValues(
    raw = NULL,
    norm = NULL,
    from = FALSE)

  ### load data
  observeEvent(input$upload, {
    mydf$raw <- rawcounts()
    mydf$norm <- norm_counts()
    mydf$from <- all(dim(rawcounts())==dim(norm_counts()))
  })
  ### clear button
  observeEvent(input$reset, {
    mydf$raw <- NULL
    mydf$norm <- NULL
    mydf$from <- FALSE
  })

  output$no_genes <- renderValueBox({
    if(mydf$from) {
      valueBox(
        value = nrow(mydf$norm),
        subtitle = "Number of Genes",
        icon = icon("dna"),
        color = "orange"
      )
    } else {
      valueBox(
        value = NA,
        subtitle = "Number of Genes",
        icon = icon("dna"),
        color = "orange"
      )
    }
  })

  output$uploadnote <- renderText({
    if(mydf$from) return("Find datasets. Begin analysis.")
    if(!mydf$from) return("No qualified datasets.")
    })


  output$no_cells <- renderValueBox({

    if(mydf$from) {
      valueBox(
        value = ncol(mydf$norm),
        subtitle = "Number of Cells",
        icon = icon("table"),
        color = "blue"
      )
    } else {
      valueBox(
        value = NA,
        subtitle = "Number of Cells",
        icon = icon("table"),
        color = "blue"
      )
    }
  })


  # Test DC in HighDim:
  step1_test2 <- reactive({
    if(!mydf$from) return(NULL)
    norm_mat <- mydf$norm
    raw_mat <- mydf$raw
    dist_mat <- parDist(t(norm_mat), method = "manhattan")
    HD_DCClusterscheck(dist_mat=dist_mat, rawcounts=raw_mat, clust.max = 5)
  })

  output$step1_dc <- renderText({
    if(!mydf$from & is.null(step1_test2())) return(NULL)
    ifelse(step1_test2()$ifConnected,
           "No. We didn't detect disconnected clusters.",
           "Yes. Upon detecting disconnected clusters, it signified the presence of diverse cell types,
           making trajectory fitting inappropriate. To aid in defining these diverse cell types,
           we present the results of the differential expression (DE) analysis between clusters.")
  })

  # DE in DC clusters:
  step1_de <- reactive({
    if(!step1_test2()$ifConnected) {
      raw_mat <- mydf$raw
      cls <- step1_test2()$Clusters
      df <- Escort::DE_seurat(rawcounts=raw_mat, cls=cls)
      return(df)
    }
  })

  output$dc_de_tb <-  DT::renderDT({
    if(!mydf$from) return(NULL)
    if(is.null(step1_test2())) return(NULL)
    if(step1_test2()$ifConnected) return(NULL)
    datatable(step1_de(),rownames = TRUE, filter = 'top')%>%
      DT::formatStyle(names(step1_de()),lineHeight='80%')
  })


  # Test Homogeneous:
  step1_test1 <- reactive({
    if(!mydf$from) return(NULL)
    norm_mat <- mydf$norm
    step1_testHomogeneous(norm_counts=norm_mat, num.sim = 1000)
  })

  output$step1_homogeneous <- renderText({
    if(!mydf$from) return(NULL)
    ifelse(step1_test1()$signal_pct>=0.46,
           "No. We could detect the trajectory signal.",
           "Yes. In the absence of a detected trajectory signal, it suggested the presence of a homogeneous dataset,
           rendering trajectory fitting inappropriate. To support this assessment, we present the results of
           highly variable genes (HVGs) analysis and GO Enrichment Analysis focusing on cell cycle.")
  })

  # HVGs in Homogeneous cells:
  step1_hvgs <- reactive({
    if(step1_test1()$signal_pct<0.46) {
      norm_mat <- mydf$norm
      df <- Escort::HVGs_quick(norm_counts=norm_mat)
      return(df)
    }
  })

  output$homo_hvgs_tb <-  DT::renderDT({
    if(!mydf$from) return(NULL)
    if(is.null(step1_test1())) return(NULL)
    if(step1_test1()$signal_pct>=0.46) return(NULL)
    datatable(step1_hvgs(),rownames = TRUE, filter = 'top')%>%
      DT::formatStyle(names(step1_de()),lineHeight='80%')
  })

  # HVGs in GO enrichment:
  step1_go <- reactive({
    if(is.null(input$go_info)) return(NULL)
    if(step1_test1()$signal_pct<0.46) {
      norm_mat <- mydf$norm
      df <- Escort::HVGs_GO(norm_counts=norm_mat, OrgDb = input$go_info)
      return(df)
    }
  })

  output$homo_go_txt <- renderText({
    if(!mydf$from) return(NULL)
    if(is.null(step1_test1())) return(NULL)
    if(step1_test1()$signal_pct>=0.46) return(NULL)
    if(is.null(step1_go())) return("There is no gene overlapping between Gene Ontology (GO) sets.")
  })

  output$homo_go_tb <- DT::renderDT({
    if(!mydf$from) return(NULL)
    if(is.null(step1_test1())) return(NULL)
    if(step1_test1()$signal_pct>=0.46) return(NULL)
    if(is.null(step1_go())) return(NULL)
    datatable(step1_go(),rownames = TRUE, filter = 'top',
              options = list(
                autoWidth = TRUE,scrollX=TRUE,
                columnDefs = list(list(width = '400px', targets = c(2)),
                                  list(visible=FALSE, targets=c(6)))))%>%
      DT::formatStyle(names(step1_de()),lineHeight='80%')
  })


  # step1_decision
  step1_res <- reactive({
    if(!mydf$from) return(NULL)
    step1_test1()$signal_pct>=0.46 && step1_test2()$ifConnecte
  })

  output$step1_decision <- renderText({
    if(!mydf$from) return(NULL)
    if(step1_res()) {
      "Go to STEP 2"
    } else {
      "Not suitable for trajectory fitting. Check results in the left column. "
    }

  })


  # Visualization
  output$step1_plot <- renderPlot({
    if(!mydf$from & is.null(step1_test1()) & is.null(step1_test2())) return(NULL)
    norm_mat <- mydf$norm
    raw_mat <- mydf$raw
    # Visualization
    K <- step1_test2()$K
    par(mfrow=c(1,2))
    library(umap)
    plotcol <- as.factor(step1_test2()$Clusters)
    dimred_umap <- umap::umap(t(norm_mat))$layout
    library(Rtsne)
    dimred_tsne <- Rtsne::Rtsne(t(norm_mat), dims = 2)$Y
    rownames(dimred_tsne) <- rownames(t(norm_mat))

    Sys.sleep(1)
    plot(dimred_umap, col = alpha(plotcol,0.7), pch=16, main="UMAP")
    legend("topright", legend=as.character(1:K), col=c(1:K), pch=16, cex = 0.5)
    plot(dimred_tsne, col = alpha(plotcol,0.7), pch=16, main="TSNE")
    legend("topleft", legend=as.character(1:K), col=c(1:K), pch=16, cex = 0.5)

  })



  # generate the obj

  step23_obj <- reactive({
    if(is.null(input$normfile)) return(NULL)
    if(!mydf$from) return(NULL)
    # select genes:
    gene.var <- quick_model_gene_var(mydf$norm)
    genes.HVGs <- rownames(gene.var)[1:input$checkgenes]
    sub_counts <- mydf$norm[genes.HVGs,]
    # DR
    dimred <- getDR_2D(sub_counts, input$checkDR)
    # Trajectory
    set.seed(123)
    cl1 <- Mclust(dimred)$classification
    ti_out <- slingshot::slingshot(data=dimred, clusterLabels=cl1)
    rawpse <- slingshot::slingPseudotime(x=ti_out, na=T)
    pse <- as.data.frame(rawpse / max(rawpse, na.rm=TRUE))
    ls_fitLine <- lapply(slingCurves(ti_out), function(x) x$s[x$ord,])
    fitLine <- do.call(rbind, lapply(ls_fitLine, function(x) {
      df_seg <- cbind(x[-nrow(x),],x[-1,])
      colnames(df_seg) <- c("x0", "y0", "x1", "y1")
      return(df_seg)
    }))
    fitLine <- as.data.frame(fitLine)

    prepTraj(dimred, PT=rawpse, fitLine=fitLine)
  })

  # plot the trajectory
  output$trajectory_plot <- renderPlot({
    if(is.null(step23_obj())) return(NULL)
    colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
    pse <- step23_obj()$pse
    pse$Ave <- rowMeans(pse, na.rm = T)

    Sys.sleep(1)
    plotcol <- colors[cut(pse$Ave, breaks=100)]
    plot(step23_obj()$Embedding, col = scales::alpha(plotcol, 0.7), pch=16, main="Estimated PT")
    segments(x0 = step23_obj()$fitLine$x0,
             y0 = step23_obj()$fitLine$y0,
             x1 = step23_obj()$fitLine$x1,
             y1 = step23_obj()$fitLine$y1, lwd = 3)
  }, height = 450, width=450)


  # download the obj
  rf2 <- reactiveValues()
  observe({
    if(!is.null(step23_obj()))
      isolate(
        eval_obj <<- step23_obj()
      )
  })

  output$downloadTraj <- downloadHandler(
    filename = function() {
      paste0(paste(input$checkDR, input$checkgenes, input$checkTraj, input$checkcls, sep="_"), ".rds")
    },
    content = function(file) {
      saveRDS(eval_obj, file = file)
    }
  )

  # load data files
  output$obj_files <- renderTable(input$objs[,1], colnames = F)

  all_files <- reactive({
    req(input$objs)
    purrr::map(input$objs$datapath, readRDS) %>%
      purrr::set_names(input$objs$name)
  })

  # step2:
  # test DC in

  step2_test1 <- reactive({
    if(is.null(all_files)) return(NULL)
    if(is.null(step1_res())) return(NULL)
    DRLvsC <- list()
    for (i in 1:length(all_files())) {
      subls <- all_files()[[i]]
      DRLvsC[[names(all_files())[i]]] <- LD_DCClusterscheck(dist_mat=dist(subls$Embedding, method = "euclidean"), DRdims=subls$Embedding, connectedCells = 1)
    }
    return(DRLvsC)
  })

  # test similarity between HD and LD
  step2_test2 <- reactive({
    if(is.null(all_files)) return(NULL)
    if(is.null(step1_res())) return(NULL)
    simi_cells <- list()
    for (i in 1:length(all_files())) {
      subls <- all_files()[[i]]
      simi_cells[[names(all_files())[i]]] <- Similaritycheck(norm_counts=mydf$norm, dimred=subls$Embedding, Cluters=step1_test2())
    }
    return(simi_cells)
  })

  #summary:
  structure_tb <- reactive({
    if(is.null(all_files)) return(NULL)
    if(is.null(step1_res())) return(NULL)
    norm_mat <- all_files()[[1]]$Normcounts
    dc_tb <- data.frame(data=names(step2_test1()), DCcheck=sapply(step2_test1(), function(x) x$ifConnected))
    simi_tb <- data.frame(data=names(step2_test2()), SimiRetain=sapply(step2_test2(), function(x) x$GoodRate))
    merge(dc_tb, simi_tb, by="data")
  })

  output$step2_structuretb <- renderTable({
    if(is.null(all_files)) return(NULL)
    if(is.null(step1_res())) return(NULL)
    df <- structure_tb()
    df$SimiRetain <- df$SimiRetain
    df$" " <- ifelse(df$DCcheck, "√", "")

    Sys.sleep(1)
    df[order(df$data), ]
  }, digits = 3)


  # test GOF
  step2_test4 <- reactive({
    if(is.null(all_files)) return(NULL)
    if(is.null(step1_res())) return(NULL)
    gof_eval <- list()
    for (i in 1:length(all_files())) {
      subls <- all_files()[[i]]
      gof_eval[[names(all_files())[i]]] <- GOFeval(subls$Embedding)
    }
    return(gof_eval)
  })

  output$step2_spreadtb <- renderTable({
    if(is.null(all_files)) return(NULL)
    if(is.null(step1_res())) return(NULL)
    gof_tb <- data.frame(data=names(step2_test4()), GOF=sapply(step2_test4(), function(x) x$occupiedRate))

    Sys.sleep(1)
    gof_tb[order(gof_tb$data), ]
  }, digits = 3)

  # step3:
  # test Ushape

  step2_test3 <- reactive({
    if(is.null(all_files)) return(NULL)
    if(is.null(step1_res())) return(NULL)
    ushap_eval <- list()
    for (i in 1:length(all_files())) {
      subls <- all_files()[[i]]
      ushap_eval[[names(all_files())[i]]] <- UshapeDetector(subls)
    }
    return(ushap_eval)
  })

  output$step3_res <- renderTable({
    if(is.null(all_files)) return(NULL)
    if(is.null(step1_res())) return(NULL)
    ushap_tb <- data.frame(data=names(step2_test3()), USHAPE=sapply(step2_test3(), function(x) x$Ambpct))

    Sys.sleep(1)
    ushap_tb[order(ushap_tb$data), ]
  }, digits = 3)

  # conclusion
  # combine all results
  final_tb <- reactive({
    if(is.null(all_files)) return(NULL)
    if(is.null(step1_res())) return(NULL)
    df <- structure_tb()
    if(!any(df$DCcheck)) return(NULL)
    # simi_tb <- data.frame(data=names(step2_test2()), SimiRetain=sapply(step2_test2(), function(x) x$GoodRate))
    ushap_tb <- data.frame(data=names(step2_test3()), USHAPE=sapply(step2_test3(), function(x) x$Ambpct))
    gof_tb <- data.frame(data=names(step2_test4()), GOF=sapply(step2_test4(), function(x) x$occupiedRate))
    df2 <- merge(ushap_tb, gof_tb, by="data")
    alldf <- merge(df, df2, by="data")
    # alldf$note[is.na(alldf$note)] <- rep(" ", sum(is.na(alldf$note)))
    rownames(alldf) <- alldf$data
    alldf[,c("DCcheck", "SimiRetain", "GOF", "USHAPE")]
  })




  res_tb <- reactive({
    if(is.null(all_files)) return(NULL)
    if(is.null(final_tb())) return(NULL)
    if(is.null(step1_res())) return(NULL)
    scoredf <- final_tb()
    final_df <- calcScore(scoredf)
    return(final_df)
  })

  output$final_res <- renderTable({
    if(is.null(all_files)) return(NULL)
    if(is.null(res_tb())) return(NULL)
    if(is.null(step1_res())) return(NULL)
    res_tb()[,c("Row.names", "DCcheck", "SimiRetain", "GOF", "USHAPE", "score","ranking", "decision", "note")]
  }, digits = 3)


  output$final_plot <- renderPlot({
    if(is.null(all_files)) return(NULL)
    if(is.null(res_tb())) return(NULL)
    if(is.null(step1_res())) return(NULL)

    df <- res_tb()
    df <- df[!is.na(df$score),]

    df <- df[order(df$score, decreasing = T), ]
    data_plt <- head(df$Row.names, 6)

    colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)

    Sys.sleep(1)
    par(mfrow = c(2, 3))
    for (i in data_plt) {
      subls <- all_files()[[i]]
      pse <- subls$pse
      pse$Ave <- rowMeans(pse, na.rm = T)
      plotcol <- colors[cut(pse$Ave, breaks=100)]
      plot(subls$Embedding, col = scales::alpha(plotcol, 0.7), pch=16, main=i)
      segments(x0 = subls$fitLine$x0,
               y0 = subls$fitLine$y0,
               x1 = subls$fitLine$x1,
               y1 = subls$fitLine$y1, lwd = 3)
    }
  })


}
