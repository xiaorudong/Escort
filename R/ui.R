library(shinydashboard)
library(shiny)
library(shinyjs)
library(shinycssloaders)


js <- '.nav-tabs-custom .nav-tabs li.active {
    border-top-color: #f39c12;
}"'


ui <- dashboardPage(
  dashboardHeader(title = "Escort"),
  dashboardSidebar(
    sidebarMenu(
      id = "tabs",
      menuItem("About", tabName = "home", icon = icon("house-user")),
      menuItem("Step 1", tabName = "step1"),
      menuItem("Generate embeddings", tabName = "generate_embeddings"),
      menuItem("Step 2", tabName = "step2"),
      menuItem("Step 3", tabName = "step3"),
      menuItem("Summary", tabName = "conclusion")
      #
      # menuItem("Start Analysis", tabName = "app", icon = icon("circle-play"),
      #          startExpanded = TRUE,
      #          # menuSubItem("Load Data", tabName = "load"),
      #          menuSubItem("Step 1", tabName = "step1"),
      #          menuSubItem("Downstream Analysis", tabName = "obj"),
      #          menuSubItem("Step 2", tabName = "step2"),
      #          menuSubItem("Step 3", tabName = "step3"),
      #          menuSubItem("Conclusion", tabName = "conclusion"))
    )
  ),
  dashboardBody(
    useShinyjs(),
    tags$style(HTML(".main-sidebar { font-size: 20px!important; }
                   .treeview-menu>li>a { font-size: 20px!important; }")),
    tags$head(tags$style(HTML("
      .main-header .logo {
        font-weight: bold;
        font-size: 24px;
      }
      .download_step_1_button {
        background-color: #0073b7;
        color: #ffffff;
        border: none;
      }
      .download_step_1_button:hover {
        background-color: #005a90;
        color: #ffffff;
      }
      .download_step_1_button:active {
        background-color: #005a90 !important;
        color: #ffffff !important;
      }
      .download_step_1_button:focus {
        background-color: #0073b7;
        color: #ffffff;
      }
    "))),
    tags$head(tags$style("
   body {
      font-size: 16px;
   }")),
    tags$style(js),
    tags$head(includeHTML("analytics.html")),
    tabItems(
      tabItem(
        tabName = "home",
        h3(strong("Welcome to Escort!")),
        fluidRow(
          column(
            12,
            "Escort is a framework that evaluates
                       various data processing decisions in terms of their effect on trajectory inference
                       with single-cell RNA-seq data. Escort guides users
                       through a trajectory analysis by providing evaluations of embeddings, which
                       represent combinations of analysis choices including feature selection,
                       dimension reduction, normalization, and/or trajectory inference-specific hyperparameters.",
            br(),
            br(),
            "In Step 1, Escort will assess the evidence of a trajectory signal in the dataset.
                       Sometimes data are not suitable for trajectory analysis,
                       for example, when cells come from biologically distinct clusters or have insufficient heterogeneity. In these cases,
                       Escort will alert the user and offers guidance to further investigate the appropriateness of trajectory analysis.",
            br(),
            br(),
            "In Step 2, Escort will compare various embeddings, specifically in terms of how well they preserve cellular relationships and the
                       distribution of cells in the embedding.",
            br(),
            br(),
            "In Step 3, Escort evaluates how well a specific trajectory inference method
                       interacts with a given embedding. This allows for evaluaton of additional graph structures used by specific methods and consideration of method-specific parameters.

                       Finally, Escort provides an overall score for each option, as well as, a classification to
                       help researchers select more optinal analysis choices for inferring a trajectory from their data.",
            br(),
            br(),
            "For any questions or issues, please submit a comment to our ", a(
              href = "https://github.com/xiaorudong/Escort/issues",
              target = "_blank",
              "GitHub issues page"
            ), ".",
            br(),
            br(),
            "Additional explanations are provided in our vignette for ", a(
              href = "https://rbacher.rc.ufl.edu/escort/docs-escort/vignettes.html",
              target = "_blank",
              "shinyEscort"
            ), "."
          )

          # column(2, imageOutput("home_img"))
        )
      ),
      tabItem(
        tabName = "step1",
        fluidRow(
          column(
            width = 6,
            box(
              width = NULL, status = "primary", title = strong("Step 1: Detecting existence of a trajectory signal"),
              "Before fitting a trajectory, Escort determines whether a given dataset is appropriate for trajectory
                            analysis. There are two common scenarios that may render trajectory analysis inappropriate or require
                            more careful consideration: ",
              br(),
              " - Datasets that contain distinct/disjoint cell types", br(),
              " - Datasets in which the cells are too homogeneous/similar",
              br(),
              br(),
              "Disjoint cell types evaluation: If this module fails, Escort will perform cluster-specific differential expression testing and display the top cluster-specific
                           genes. Users should further examine this list to determine whether fitting a trajectory that connects these particular cell types is biologically reasonable.
                           If so, then users should re-examine whether sufficient intermediate cells exist between these two groups, or whether batch effects have been exist and have been addressed.",
              br(),
              br(),
              "Homogenous cells evaluation: If this module fails, Escort will report the top highly variable genes in the dataset along with their enrichments.
                           Users should further examine this list to re-examine whether an underlying trajectory is appropriate. If so, then this output should be used to identify what
                           other biological processes may be overriding or diluting the biological signal of interest (e.g. cell cycle) or excessive signal from ribosomal or mitochondrial genes."
            ),
            tabBox(
              title = "", width = NULL, id = "dc",
              tabPanel("About", strong("Disjoint cell types detected?"), textOutput(outputId = "step1_dc"), hr(), downloadButton("downloadStep1de", "Download DE", class = "download_step_1_button"), ),
              tabPanel("DE", DT::DTOutput(outputId = "dc_de_tb") %>% withSpinner(color = "#FAD02C"))
            ),
            tabBox(
              title = "", width = NULL, id = "homo",
              tabPanel("About", strong("Homogeneous cells detected?"), textOutput(outputId = "step1_homogeneous"), hr(), downloadButton("downloadStep1hvg", "Download HVG", class = "download_step_1_button"), downloadButton("downloadStep1go", "Download GO", class = "download_step_1_button")),
              tabPanel("HVGs", DT::DTOutput(outputId = "homo_hvgs_tb") %>% withSpinner(color = "#FAD02C")),
              tabPanel(
                "GO",
                helpText("Please select a database that contains information about the genome and gene annotations of a particular organism."),
                radioButtons("go_info", "Data From",
                  choices = list("Human" = "org.Hs.eg.db", "Mouse" = "org.Mm.eg.db"), selected = 1
                ),
                textOutput(outputId = "homo_go_txt"),
                DT::DTOutput(outputId = "homo_go_tb") %>% withSpinner(color = "#FAD02C")
              )
            )
          ),
          column(
            width = 6,
            h4(strong("Upload scRNA-seq datasets:")),
            "Please upload either .csv files or .rds files for both raw data and normalized data. More information on how to prepare these files from a Seurat object or SingleCellExperiment is provided",
            a(
              href = "https://rbacher.rc.ufl.edu/escort/docs-escort/vignettes.html",
              target = "_blank",
              "here."
            ),
            fluidRow(
              column(
                width = 6,
                br(),
                actionButton("upload_example_data", label = "Load Example Data"),
                br(),
                br(),
                # read raw data after QC
                fileInput("rawfile", label = "Raw Data", buttonLabel = "File", accept = c(".csv", ".rds")),
                # read norm data after QC
                fileInput("normfile", label = "Normalized Data", buttonLabel = "File", accept = c(".csv", "rds")),
                actionButton("upload", "Upload"),
                actionButton("reset", "Clear"), # Clear input dataset#
                br(),
                br(),
                textOutput(outputId = "uploadnote"),
              ),
              column(
                width = 6,
                br(),
                # h4(strong("Data Info:")),
                valueBoxOutput("no_cells", width = NULL),
                valueBoxOutput("no_genes", width = NULL)
              )
            ),
            fluidRow(
              column(
                width = 12,
                br(),
                box(
                  width = NULL, status = "warning",
                  span(uiOutput(outputId = "step1_decision"), style = "font-size:18px; font-style:bold; text-align:left"),
                  hr(),
                  downloadButton("downloadStep1Results", "Download Step 1 Results", class = "download_step_1_button"),
                  hr(),
                  actionButton("go_to_generate_embeddings", label = "Go to next step"),
                )
              )
            )
          )
        ),
      ),
      tabItem(
        tabName = "generate_embeddings",
        fluidRow(
          column(
            3,
            h4(strong("Generate embeddings:")),
            "Embeddings can be generated invidivually here (in combination with the preferred trajectory method used in Step 3), or users can generate
                      their own embeddings following the workflow in the ", a(
              href = "https://rbacher.rc.ufl.edu/escort/docs-escort/generate_embedding_objects_vignette.html",
              target = "_blank",
              "Generating data objects for Escort R/Shiny Vignette"
            ), ".",
            numericInput("checkgenes", "The number of selected HVGs", value = 100, min = 0),
            # choose DR method
            checkboxGroupInput(
              "drMethods",
              "Dimension reduction method",
              choices = c(
                "MDS" = "MDS",
                "TSNE" = "TSNE",
                "UMAP" = "UMAP"
              ),
              selected = c("MDS", "TSNE", "UMAP")
            ),
            h4(strong("Fit preliminary trajectory:")),
            # choose Trajectory methods
            selectInput("checkTraj", "Trajectory method", choices = c("Slingshot" = "Slingshot")),
            br(),
            downloadButton(outputId = "downloadTraj", label = "Download .rds object")
          ),
          column(3, plotOutput("trajectory_plot_MDS", height = "auto") %>% withSpinner(color = "#FAD02C", id = "trajectory_plot_MDS_spinner")),
          column(3, plotOutput("trajectory_plot_TSNE", height = "auto") %>% withSpinner(color = "#FAD02C", id = "trajectory_plot_TSNE_spinner")),
          column(3, plotOutput("trajectory_plot_UMAP", height = "auto") %>% withSpinner(color = "#FAD02C", id = "trajectory_plot_UMAP_spinner")),
          br(),
          br(),
          column(
            9, # New column for the upload button
            br(),
            br(),
            fileInput("normalfile", label = "Optional: Upload the Normalized Data", buttonLabel = "Upload", accept = c(".csv", "rds"))
          )
        )
      ),
      tabItem(
        tabName = "step2",
        fluidRow(
          column(
            width = 9,
            box(
              width = NULL, status = "primary", title = strong("Step 2: Evaluating the trajectory characteristics of embeddings"),
              "Next, Escort identifies preferred embeddings for performing trajectory inference.
                       We will evaluate three characteristics of each embedding: ", br(),
              " - The retention of inter-cellular relationships.", br(),
              " - The preservation of similarity relationships.", br(),
              " - Distribution of cells in the embedding space."
            ),
            box(
              width = NULL, title = "Inter-cellular relationships", status = "warning",
              tableOutput(outputId = "step2_celltb") %>% withSpinner(color = "#FAD02C")
            ),
            box(
              width = NULL, title = "Preservation of similarity relationships", status = "warning",
              tableOutput(outputId = "step2_simitb") %>% withSpinner(color = "#FAD02C")
            ),
            box(
              width = NULL, title = "Cell spread ", status = "warning",
              tableOutput(outputId = "step2_spreadtb") %>% withSpinner(color = "#FAD02C")
            )
          ),
          column(
            width = 3,
            h4(strong("Load all embeddings: (multiple allowed)")),
            fileInput("objs", label = NULL, buttonLabel = "Upload .rds file", accept = c(".rds"), multiple = TRUE),
            tableOutput("obj_files"),
            fileInput("step1ResultsUpload", buttonLabel = "Upload .rds file", "Upload Step 1 Results", accept = ".rds")
          )
        )
      ),
      tabItem(
        tabName = "step3",
        fluidRow(
          box(
            width = 12, status = "primary", title = strong("Step 3: Quantifying trajectory fitting performance"),
            "Embeddings are also evaluated in the context of a trajectory inference method.
                A preliminary trajectory in inferred and Escort assesses the proportion of cells having an ambiguous projection
                to the trajectory. For example, trajectories in a U-shape tend to
                be inaccurate because some cells will have map to opposing pseudotimes with equal probability."
          )
        ),
        fluidRow(
          box(
            width = 12, title = "Percentage of ambiguous cells", status = "warning",
            tableOutput(outputId = "step3_res") %>% withSpinner(color = "#FAD02C")
          )
        )
      ),
      tabItem(
        tabName = "conclusion",
        fluidRow(
          box(
            width = 12, status = "primary", title = strong("Escort suggestions for embedding selection"),
            "Below is a table with each embedding's rating according to their overall performance. Embeddings with a score larger than zero are
                  Recommended for trajectoy inference.",
            hr(),
            tableOutput(outputId = "final_res") %>% withSpinner(color = "#FAD02C")
          )
        ),
        downloadButton(outputId = "table", label = "Download the table"),
        hr(),
        fluidRow(
          box(
            width = 12, status = "warning",
            plotOutput("final_plot", height = "650px", width = "900px") %>% withSpinner(color = "#FAD02C")
          )
        )
      )
    )
  )
)
