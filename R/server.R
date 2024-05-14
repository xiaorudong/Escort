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
    tryCatch({
      if (grepl("\\.csv$", tolower(input$rawfile$name), ignore.case = TRUE)) {
        mat <- read.csv(input$rawfile$datapath, row.names = 1)
        as.matrix(mat)
      } else if (grepl("\\.rds$", tolower(input$rawfile$name), ignore.case = TRUE)) {
        mat <- readRDS(input$rawfile$datapath)
        if(!is.matrix(mat)) {
          stop("The RDS file must contain a matrix.")
        }
        mat
      } else {
        stop("Unsupported file format. Please upload a .csv or .rds file.")
      }
    }, error = function(e) {
      showNotification(e$message, type = "error")
      NULL
    })
  })
  
  

  normcounts <- reactive({
    req(input$normfile)
    tryCatch({
      if (grepl("\\.csv$", tolower(input$normfile$name), ignore.case = TRUE)) {
        mat <- read.csv(input$normfile$datapath, row.names = 1)
        as.matrix(mat)
      } else if (grepl("\\.rds$", tolower(input$normfile$name), ignore.case = TRUE)) {
        mat <- readRDS(input$normfile$datapath)
        if(!is.matrix(mat)) {
          stop("The RDS file must contain a matrix.")
        }
        mat
      } else {
        stop("Unsupported file format. Please upload a .csv or .rds file.")
      }
    }, error = function(e) {
      showNotification(e$message, type = "error")
      NULL
    })
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

  ### load example data 
  observeEvent (input$upload_example_data,{
    isolate({
      rawcount_linear_splatter <- read.csv("data/rawcount_linear_splatter.csv", row.names = 1)
      normcount_linear_splatter <- read.csv("data/normcount_linear_splatter.csv", row.names = 1)

      as.matrix(rawcount_linear_splatter)
      as.matrix(normcount_linear_splatter)

      mydf$raw <- rawcount_linear_splatter
      mydf$norm <- normcount_linear_splatter
      mydf$from <- all(dim(normcount_linear_splatter)==dim(normcount_linear_splatter))
      
			showNotification("Example data is loading. Please wait.", type="default", duration=8)
		})
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
    Escort::HD_DCClusterscheck(normcounts=norm_mat, rawcounts=raw_mat, clust.max = 5)
  })
  
 

  output$step1_dc <- renderText({
    if(!mydf$from & is.null(step1_test2())) return(NULL)
    ifelse(step1_test2()$ifConnected,
           "No. We didn't detect disconnected clusters.",
           "Yes. Upon detecting disconnected clusters, it signified the presence of diverse cell types,
           making trajectory fitting inappropriate. To aid in defining these diverse cell types,
           we present the results of the differential expression (DE) analysis between clusters.")
  })
  
  #save the result from step 1:
  output$downloadStep1Results <- downloadHandler(
    # Define the filename of the download
    filename = function() {
      paste("step_1_results", Sys.Date(), ".rds", sep = "")
    },
    content = function(file) {
      req(step1_test2())
      saveRDS(step1_test2(), file)
    }
  )
  

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
    DT::datatable(step1_de(),rownames = F, filter = 'top')%>%
      DT::formatStyle(names(step1_de()),lineHeight='80%')
  })


  # Test Homogeneous:
  step1_test1 <- reactive({
    if(!mydf$from) return(NULL)
    norm_mat <- mydf$norm
    Escort::step1_testHomogeneous(normcounts=norm_mat, num.sim = 1000)
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
    DT::datatable(step1_hvgs(),rownames = TRUE, filter = 'top')%>%
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
    DT::datatable(step1_go(),rownames = TRUE, filter = 'top',
              options = list(
                autoWidth = TRUE,scrollX=TRUE,
                columnDefs = list(list(width = '400px', targets = c(2)),
                                  list(visible=FALSE, targets=c(6)))))%>%
      DT::formatStyle(names(step1_de()),lineHeight='80%')
  })


  # step1_decision
  step1_res <- reactive({
    if(!mydf$from)return(NULL)
    step1_test1()$signal_pct>=0.46 && (step1_test2()$ifConnected)
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
    dimred_umap <- Escort::getDR_2D(norm_mat, "UMAP")
    # dimred_umap <- umap::umap(t(norm_mat))$layout
    # library(Rtsne)
    dimred_tsne <- Escort::getDR_2D(norm_mat, "TSNE")
    # dimred_tsne <- Rtsne::Rtsne(t(norm_mat), dims = 2)$Y
    # rownames(dimred_tsne) <- rownames(t(norm_mat))

    Sys.sleep(1)
    plot(dimred_umap, col = alpha(plotcol,0.7), pch=16, main="UMAP")
    legend("topright", legend=as.character(1:K), col=c(1:K), pch=16, cex = 0.5)
    plot(dimred_tsne, col = alpha(plotcol,0.7), pch=16, main="TSNE")
    legend("topleft", legend=as.character(1:K), col=c(1:K), pch=16, cex = 0.5)

  })
  
  
  #optional for uploading norm data in embedding
  normalcounts <- reactive({
    req(input$normalfile)
    tryCatch({
      if (grepl("\\.csv$", tolower(input$normalfile$name), ignore.case = TRUE)) {
        mat <- read.csv(input$normalfile$datapath, row.names = 1)
        as.matrix(mat)
      } else if (grepl("\\.rds$", tolower(input$normalfile$name), ignore.case = TRUE)) {
        mat <- readRDS(input$normalfile$datapath)
        if(!is.matrix(mat)) {
          stop("The RDS file must contain a matrix.")
        }
        mat
      } else {
        stop("Unsupported file format. Please upload a .csv or .rds file.")
      }
    }, error = function(e) {
      showNotification(e$message, type = "error")
      NULL
    })
  })
  
  

  observeEvent(input$normalfile, {
    optional_normal_file <- normalcounts()
    if(!is.null(optional_normal_file)) {
      mydf$norml <- optional_normal_file
      mydf$readyForEmbedding <- TRUE  
    } else {
      mydf$readyForEmbedding <- FALSE
    }
  })
  
  #optional safegenes
  safegenes <- reactive({
    
    req(isTRUE(mydf$readyForEmbedding) || !is.null(mydf$norm))
    
    
    currentData <- if (!is.null(mydf$norm)) {
      mydf$norm
    } else if (!is.null(mydf$norml)) {
      mydf$norml
    } else {
      return(NULL)  
    }
    
    
    max_genes <- nrow(currentData)
    
    validated_genes <- NULL
    
    if (is.na(input$checkgenes) || is.null(input$checkgenes)) {
      validated_genes <- 10
      showNotification("No inpit detected or input is invalid. Defaulting to 10 genes.", type = "message")
    } else {
      validated_genes <- min(max(input$checkgenes, 10), max_genes)
    
    # Notifications based on the outcome of `validated_genes`
    if (input$checkgenes != validated_genes) {
      if (input$checkgenes < 10) {
        showNotification("Input is below the minimum allowed number (10). Adjusted to 10.", type = "message")
      } else if (input$checkgenes > max_genes) {
        showNotification("Input exceeds the maximum number of genes. Adjusted to maximum available.", type = "message")
      }
    }
    }
     return(validated_genes)
  })
  
  

  
  step23_obj <- reactive({
    req(isTRUE(mydf$readyForEmbedding) || !is.null(mydf$norm))
    currentData <- if (!is.null(mydf$norm)) mydf$norm else mydf$norml
    tryCatch({
      # Select genes
      gene.var <- Escort::quick_model_gene_var(currentData)
      genes.HVGs <- rownames(gene.var)[1:safegenes()]
      sub_counts <- currentData[genes.HVGs,]

      #all types of dr method embeddings to be entered in the list
      dr_methods <- c("MDS", "TSNE", "UMAP")
      #This list will contain list of objects (objects returned by the porpTraj method)
      #This list will also be returned
      list_object_to_return <- list()
      #Loop through all the selected DR methods
      for(dr_method in dr_methods) {
        dimred <- Escort::getDR_2D(sub_counts, dr_method)
        # Trajectory
        set.seed(123)
        cl1 <- mclust::Mclust(dimred)$classification
        ti_out <- slingshot::slingshot(data=dimred, clusterLabels=cl1)
        rawpse <- slingshot::slingPseudotime(x=ti_out, na=T)
        pse <- as.data.frame(rawpse / max(rawpse, na.rm=TRUE))
        ls_fitLine <- lapply(slingshot::slingCurves(ti_out), function(x) x$s[x$ord,])
        fitLine <- do.call(rbind, lapply(ls_fitLine, function(x) {
          df_seg <- cbind(x[-nrow(x),],x[-1,])
          colnames(df_seg) <- c("x0", "y0", "x1", "y1")
          return(df_seg)
        }))
        fitLine <- as.data.frame(fitLine)
        returned_object <- Escort::prepTraj(dimred, PT=rawpse, fitLine=fitLine)
        list_object_to_return[[dr_method]] <- returned_object
      }
      return(list_object_to_return)
    }, error = function(e) {
      showNotification(paste("The input number is invalid, minimum is 10:", e$message), type = "error")
      return(NULL)
    })
  })
  
  
  # colors obtained from brewer.pal(11,'Spectral')[-6]
  brewCOLS <- c("#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#E6F598", "#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2")
  
  #utility function to generate plot for embeddings
  generate_embedding_plot <- function(embedding_name, embedding_object) {
    if (is.null(embedding_object)) return(NULL)
    
    colors <- colorRampPalette(brewCOLS)(100)
    pse <- embedding_object$pse
    pse$Ave <- rowMeans(pse, na.rm = TRUE)
    
    Sys.sleep(1)
    plotcol <- colors[cut(pse$Ave, breaks = 100)]
    plot(embedding_object$Embedding, col = scales::alpha(plotcol, 0.7), pch = 16, main = paste("Estimated PT -", embedding_name))
    segments(x0 = embedding_object$fitLine$x0,
            y0 = embedding_object$fitLine$y0,
            x1 = embedding_object$fitLine$x1,
            y1 = embedding_object$fitLine$y1, lwd = 3)
  }

  #Plot trajectories
  output$trajectory_plot_MDS <- renderPlot({
    temp_object <- step23_obj()
    # as the step23_obj() now contains MDS, TSNE and UMAP embeddings object
    # we pass the object which we need to plot for in the 
    # generate_embedding_plot utility function
    generate_embedding_plot("MDS", temp_object$MDS)
  }, height = 300, width = 300)

  output$trajectory_plot_TSNE <- renderPlot({
    temp_object <- step23_obj()
    # as the step23_obj() now contains MDS, TSNE and UMAP embeddings object
    # we pass the object which we need to plot for in the 
    # generate_embedding_plot utility function
    generate_embedding_plot("TSNE", temp_object$TSNE)
  }, height = 300, width = 300)

  output$trajectory_plot_UMAP <- renderPlot({
    temp_object <- step23_obj()
    # as the step23_obj() now contains MDS, TSNE and UMAP embeddings object
    # we pass the object which we need to plot for in the 
    # generate_embedding_plot utility function
    generate_embedding_plot("UMAP", temp_object$UMAP)
  }, height = 300, width = 300)

  
  output$downloadTraj <- downloadHandler(
    filename = function() {
      paste0(paste(safegenes(), input$checkTraj, sep="_"), ".rds")
    },
    content = function(file) {
      all_dr_types_embeddings_list <- step23_obj()
      
      #get all the user selected DR methods and store them in a list
      dr_methods_list <- c()
      if (!is.null(input$drMethods)) {
        dr_methods_list <- strsplit(input$drMethods, ", ")
      }
      to_download <- list()

      #loop through the umbrella list all_dr_types_embeddings_list
      #which contains embedding object for all 3 types
      #and only filter out the ones selected by the user
      #push those items to a list to_download which will finally be downloaded
      if(!is.null(all_dr_types_embeddings_list))
        for(dr_method in dr_methods_list) {
          to_download[[dr_method]] <- all_dr_types_embeddings_list[[dr_method]]
        }
      saveRDS(to_download, file = file)
    }
  )

  observe({
      #by default initially keep the download button disabled
      shinyjs::disable("downloadTraj")
      
      dr_methods_string <- input$drMethods
      is_dr_methods_checkbox_empty <- is.null(dr_methods_string)
      #check if input methods are selected and show a warning based on that
      if(is_dr_methods_checkbox_empty) {
        showNotification("At least one DR method needs to be selected.", type = "error", duration = 10)
        shinyjs::disable("downloadTraj")
        return(NULL)
      }
      #if the data is still not imported/ready then download button
      #should be disabled
      if(isFALSE(mydf$readyForEmbedding) || is.null(mydf$norm)) {
        shinyjs::disable("downloadTraj")
      } else {
        # else enable the download button
        shinyjs::enable("downloadTraj")
      }
    })

  # load data files
  output$obj_files <- renderTable(input$objs[,1], colnames = F)

  all_files <- reactive({
    req(input$objs)
    list_to_return <- list()
    rds_list <- purrr::map(input$objs$datapath, readRDS) %>%
      purrr::set_names(input$objs$name)
    
    #Once we have all the .rds files uploaded, we need to split the objects
    #inside those files
    #e.g 100_Slingshot has MDS, TSNE, UMAP objects in it
    #e.g 10_Slingshot has only MDS in it
    #so we need to return a total of 4 items to the all_files reactive property
    #i.e MDS_100_Slingshot, TSNS_100_Slingshot, UMAP_100_Slingshot and MDS_10_Slingshot
    for(rds_object_name in names(rds_list)) {
      for(dr_name in names(rds_list[[rds_object_name]])) {
        updated_keyword_name <- paste(dr_name, "_", rds_object_name, sep = "")
        list_to_return[[updated_keyword_name]] <- rds_list[[rds_object_name]][[dr_name]]
      }
    }
    return(list_to_return)
  })

  HDDCcheck_upload <- reactive({
    req(input$step1ResultsUpload)
    tryCatch({
      if (grepl("\\.rds$", tolower(input$step1ResultsUpload$name), ignore.case = TRUE)) {
        mat <- readRDS(input$step1ResultsUpload$datapath)
        
        if (!is.null(mat)) {
          return(mat)
        } else {
          showNotification("The RDS file does not contain valid data.", type = "error")
          return(NULL)
        }
      } else {
        stop("Unsupported file format. Please upload a .rds file.")
      }
    }, error = function(e) {
      showNotification(e$message, type = "error")
      NULL
    })
  })
  
  
  
  
  #optional 
  observeEvent(input$step1ResultsUpload, {
    HDDCcheck_result <- HDDCcheck_upload()
    if(!is.null(HDDCcheck_result)) {
      mydf$result <- HDDCcheck_result
  }
    })

  
  

  # step2:
  # test DC in
  
  step2_test1 <- reactive({
    if(is.null(all_files)) return(NULL)
    if(is.null(step1_res()) && is.null(mydf$result)) return(NULL)
    DRLvsC <- list()
    for (i in 1:length(all_files())) {
      subls <- all_files()[[i]]
      DRLvsC[[names(all_files())[i]]] <- Escort::LD_DCClusterscheck(embedding=subls$Embedding, connectedCells = 1)
    }
    names(DRLvsC) <- names(all_files())
    return(DRLvsC)
  })
  
  # test similarity between HD and LD
  step2_test2 <- reactive({
    if(is.null(all_files)) return(NULL)
    if(is.null(step1_res()) && is.null(mydf$result)) return(NULL)
    clusters_data <- if (!is.null(mydf$result)) {
      mydf$result
    } else {
      step1_test2()
    }
    simi_cells <- list()
    for (i in 1:length(all_files())) {
      subls <- all_files()[[i]]
      simi_cells[[names(all_files())[i]]] <- Similaritycheck(dimred=subls$Embedding, clusters=clusters_data)
    }
    names(simi_cells) <- names(all_files())
    return(simi_cells)
  })
  
  #summary:
  structure_tb <- reactive({
    if(is.null(all_files)) return(NULL)
    if(is.null(step1_res()) && is.null(mydf$result)) return(NULL)
    # norm_mat <- all_files()[[1]]$Normcounts
    dc_tb <- data.frame(data=names(step2_test1()), DCcheck=sapply(step2_test1(), function(x) x$ifConnected))
    simi_tb <- data.frame(data=names(step2_test2()), SimiRetain=sapply(step2_test2(), function(x) x$GoodRate))
    merge(dc_tb, simi_tb, by="data")
  })
  
  output$step2_celltb <- renderTable({
    if(is.null(all_files)) return(NULL)
    if(is.null(step1_res()) && is.null(mydf$result)) return(NULL)
    df <- structure_tb()[,c("data", "DCcheck")]
    df$" " <- ifelse(df$DCcheck, "√", "X")
    
    Sys.sleep(1)
    df[order(df$data), ]
  }, digits = 3)
  
  output$step2_simitb <- renderTable({
    if(is.null(all_files)) return(NULL)
    if(is.null(step1_res()) && is.null(mydf$result)) return(NULL)
    df <- structure_tb()[,c("data", "SimiRetain")]
    
    Sys.sleep(1)
    df[order(df$data), ]
  }, digits = 3)
  
  
  # test GOF
  step2_test4 <- reactive({
    if(is.null(all_files)) return(NULL)
    if(is.null(step1_res()) && is.null(mydf$result)) return(NULL)
    gof_eval <- list()
    for (i in 1:length(all_files())) {
      subls <- all_files()[[i]]
      gof_eval[[names(all_files())[i]]] <- Escort::GOFeval(subls$Embedding)
    }
    return(gof_eval)
  })
  
  output$step2_spreadtb <- renderTable({
    if(is.null(all_files)) return(NULL)
    if(is.null(step1_res()) && is.null(mydf$result)) return(NULL)
    gof_tb <- data.frame(data=names(step2_test4()), GOF=sapply(step2_test4(), function(x) x$occupiedRate))
    
    Sys.sleep(1)
    gof_tb[order(gof_tb$data), ]
  }, digits = 3)
  
  
  
  # step3:
  # test Ushape
  
  step2_test3 <- reactive({
    if(is.null(all_files)) return(NULL)
    if(is.null(step1_res()) && is.null(mydf$result)) return(NULL)
    ushap_eval <- list()
    for (i in 1:length(all_files())) {
      subls <- all_files()[[i]]
      ushap_eval[[names(all_files())[i]]] <- Escort::UshapeDetector(subls)
    }
    return(ushap_eval)
  })
  
  output$step3_res <- renderTable({
    if(is.null(all_files)) return(NULL)
    if(is.null(step1_res()) && is.null(mydf$result)) return(NULL)
    ushap_tb <- data.frame(data=names(step2_test3()), USHAPE=sapply(step2_test3(), function(x) x$Ambpct))
    
    Sys.sleep(1)
    ushap_tb[order(ushap_tb$data), ]
  }, digits = 3)
  
  # conclusion
  # combine all results
  final_tb <- reactive({
    if(is.null(all_files)) return(NULL)
    if(is.null(step1_res()) && is.null(mydf$result)) return(NULL)
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
    if(is.null(step1_res()) && is.null(mydf$result)) return(NULL)
    scoredf <- final_tb()
    final_df <- Escort::calcScore(scoredf)
    return(final_df)
  })
  
  output$final_res <- renderTable({
    if(is.null(all_files)) return(NULL)
    if(is.null(res_tb())) return(NULL)
    if(is.null(step1_res()) && is.null(mydf$result)) return(NULL)
    res_tb()[,c("Row.names", "DCcheck", "SimiRetain", "GOF", "USHAPE", "score","ranking", "decision", "note")]
  }, digits = 3)
  
  
  output$table <- downloadHandler(
      filename = function() {
        paste("final_res", ".csv", sep="")
      },
      content = function(file) {
        req(res_tb()) 
        selected_data <- res_tb()[,c("Row.names", "DCcheck", "SimiRetain", "GOF", "USHAPE", "score", "ranking", "decision", "note")]
        write.csv(selected_data, file, row.names = FALSE)
      }
    )
  
  output$final_plot <- renderPlot({
    if(is.null(all_files)) return(NULL)
    if(is.null(res_tb())) return(NULL)
    if(is.null(step1_res()) && is.null(mydf$result)) return(NULL)
    
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