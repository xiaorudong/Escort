library(mclust)
library(scales)
library(slingshot)
library(DT)
library(shinyjs)
library(dplyr)
library(clusterProfiler)
library(edgeR)


server <- function(input, output) {
  # # add image:
  # output$home_img <- renderImage({
  #
  #   list(src = "images/Workflow.png",
  #        width = 600,
  #        height = 500)
  # }, deleteFile = F)

  # step0:
  rawcounts <- reactive({
    req(input$rawfile)
    mat <- read.csv(input$rawfile$datapath, row.names=1)
    as.matrix(mat)
  })

  normcounts <- reactive({
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
    mydf$norm <- normcounts()
    mydf$from <- all(dim(rawcounts())==dim(normcounts()))
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
    if(mydf$from) return("Datasets verified. Begin analysis.")
    if(!mydf$from) return("No datasets detected.")
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
    HD_DCClusterscheck(normcounts=norm_mat, rawcounts=raw_mat, clust.max = 5)
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
      df <- Escort::DE_edgeR(rawcounts=raw_mat, cls=cls)
      return(df)
    }
  })

  output$dc_de_tb <-  DT::renderDT({
    if(!mydf$from) return(NULL)
    if(is.null(step1_test2())) return(NULL)
    if(step1_test2()$ifConnected) return(NULL)
    datatable(step1_de(),rownames = F, filter = 'top')%>%
      DT::formatStyle(names(step1_de()),lineHeight='80%')
  })


  # Test Homogeneous:
  step1_test1 <- reactive({
    if(!mydf$from) return(NULL)
    norm_mat <- mydf$norm
    step1_testHomogeneous(normcounts=norm_mat, num.sim = 1000)
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
      df <- Escort::HVGs_quick(normcounts=norm_mat)
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
      df <- Escort::HVGs_GO(normcounts=norm_mat, OrgDb = input$go_info)
      return(df)
    }
  })

  output$homo_go_txt <- renderText({
    if(!mydf$from) return(NULL)
    if(is.null(step1_test1())) return(NULL)
    if(step1_test1()$signal_pct>=0.46) return(NULL)
    if(is.null(step1_go())) return("There are no genes overlapping Gene Ontology (GO) sets.")
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
    # library(umap)
    plotcol <- as.factor(step1_test2()$Clusters)
    dimred_umap <- getDR_2D(norm_mat, "UMAP")
    # dimred_umap <- umap::umap(t(norm_mat))$layout
    # library(Rtsne)
    dimred_tsne <- getDR_2D(norm_mat, "TSNE")
    # dimred_tsne <- Rtsne::Rtsne(t(norm_mat), dims = 2)$Y
    # rownames(dimred_tsne) <- rownames(t(norm_mat))

    Sys.sleep(1)
    plot(dimred_umap, col = alpha(plotcol,0.7), pch=16, main="UMAP")
    legend("topright", legend=as.character(1:K), col=c(1:K), pch=16, cex = 0.5)
    plot(dimred_tsne, col = alpha(plotcol,0.7), pch=16, main="TSNE")
    legend("topleft", legend=as.character(1:K), col=c(1:K), pch=16, cex = 0.5)

  })



  # generate the obj

  safegenes <- reactive({
    if(is.null(input$normfile)) return(NULL)
    if(!mydf$from) return(NULL)
    # select genes:
    a <- ifelse(input$checkgenes>nrow(mydf$norm), nrow(mydf$norm), input$checkgenes)
    a <- ifelse(input$checkgenes<10, 10, input$checkgenes)
    a <- ifelse(is.na(input$checkgenes), 10, input$checkgenes)
    return(a)
  })

  step23_obj <- reactive({
    if(is.null(input$normfile)) return(NULL)
    if(!mydf$from) return(NULL)
    # select genes:
    gene.var <- quick_model_gene_var(mydf$norm)
    genes.HVGs <- rownames(gene.var)[1:safegenes()]
    sub_counts <- mydf$norm[genes.HVGs,]
    # DR
    dimred <- getDR_2D(sub_counts, input$checkDR)
    # Trajectory
    set.seed(123)
    cl1 <- mclust::Mclust(dimred)$classification
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

  # colors obtained from brewer.pal(11,'Spectral')[-6]
  brewCOLS <- c("#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#E6F598", "#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2")

  # plot the trajectory
  output$trajectory_plot <- renderPlot({
    if(is.null(step23_obj())) return(NULL)
    colors <- colorRampPalette(brewCOLS)(100)
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
      paste0(paste(input$checkDR, safegenes(), input$checkTraj, sep="_"), ".rds")
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
      DRLvsC[[names(all_files())[i]]] <- LD_DCClusterscheck(embedding=subls$Embedding, connectedCells = 1)
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
      simi_cells[[names(all_files())[i]]] <- Similaritycheck(normcounts=mydf$norm, dimred=subls$Embedding, Cluters=step1_test2())
    }
    return(simi_cells)
  })

  #summary:
  structure_tb <- reactive({
    if(is.null(all_files)) return(NULL)
    if(is.null(step1_res())) return(NULL)
    # norm_mat <- all_files()[[1]]$Normcounts
    dc_tb <- data.frame(data=names(step2_test1()), DCcheck=sapply(step2_test1(), function(x) x$ifConnected))
    simi_tb <- data.frame(data=names(step2_test2()), SimiRetain=sapply(step2_test2(), function(x) x$GoodRate))
    merge(dc_tb, simi_tb, by="data")
  })

  output$step2_celltb <- renderTable({
    if(is.null(all_files)) return(NULL)
    if(is.null(step1_res())) return(NULL)
    df <- structure_tb()[,c("data", "DCcheck")]
    df$" " <- ifelse(df$DCcheck, "√", "X")

    Sys.sleep(1)
    df[order(df$data), ]
  }, digits = 3)

  output$step2_simitb <- renderTable({
    if(is.null(all_files)) return(NULL)
    if(is.null(step1_res())) return(NULL)
    df <- structure_tb()[,c("data", "SimiRetain")]

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

    # colors obtained from brewer.pal(11,'Spectral')[-6]
    brewCOLS <- c("#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#E6F598", "#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2")

    colors <- colorRampPalette(brewCOLS)(100)

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
