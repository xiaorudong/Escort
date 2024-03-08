
library(shinydashboard)
library(shiny)
library(shinyjs)
library(shinycssloaders)
library(magrittr)

js <- '.nav-tabs-custom .nav-tabs li.active {
    border-top-color: #f39c12;
}"'


ui <- dashboardPage(
  dashboardHeader(title = "Escort"),

  dashboardSidebar(
    sidebarMenu(id="tabs",
                menuItem("About", tabName = "home", icon = icon("house-user")),
                menuItem("Step 1", tabName = "step1"),
                menuItem("Generate embeddings", tabName = "obj"),
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

    tags$head(tags$style(HTML('
      .main-header .logo {
        font-weight: bold;
        font-size: 24px;
      }
    '))),

    tags$head(tags$style('
   body {
      font-size: 16px;
   }'
    )),

    tags$style(js),

    tabItems(
      tabItem(tabName = "home",
              h3(strong("Welcome to Escort!")),
              fluidRow(
                column(8,
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
                       "For any questions or issues, please submit a comment to our ", a(href = "https://github.com/xiaorudong/Escort/issues",
                         "GitHub issues page"), ".",
                       br(),
                       br(),
                       "Additional explanations are provided in our vignette for ", a(href = "https://www.rhondabacher.com/docs-escort/shiny_vignette.html",
                         "shinyEscort"), ".")

                # column(2, imageOutput("home_img"))
                )),

      tabItem(tabName = "step1",

              fluidRow(
                column(width = 6,
                       box(width=NULL, status = "primary", title = strong("Step 1: Detecting trajectory existence"),
                           "In the first step of Escort, evidence of a trajectory signal is assessed in two scenarios: ",
                           br(),
                           " - Dataset contains distinct cell types",  br(),
                           " - Cells are too homogeneous",
                           br(),
                           br(),
                           "If the distinct cell type module fails, additional information about cluster-specific differentially expressed genes is provided. Users should examine
                           whether fitting a trajectory that connects these cell types is biologically reasonable. If so, then users should re-examine whether intermediate cell types
                           exist, presence of batch effects, and choice of normalization method.",
                           br(),
                           br(),
                           "If the homogenous cells module fails, the top highly variable genes are shown along with their enrichments. Again, the appropriateness of an underlying trajectory should
                           be considered. To proceed, users should investigate whether other processes could be overriding the biological signal of interest (e.g. cell cycle) or excessive
                            signal from ribosomal or mitochondrial genes.")),
                  column(width=6,
                         h4(strong("Upload scRNA-seq datasets:")),
                         "Note: please upload .csv files or .rds files for raw data and normalized data.",
                         fluidRow(
                           column(width = 6,
                                br(),
                                # read raw data after QC
                                fileInput("rawfile", label = "Raw Data", buttonLabel = "Upload", accept = c(".csv")),
                                # read norm data after QC
                                fileInput("normfile", label = "Norm Data", buttonLabel = "Upload", accept = c(".csv")),
                                actionButton("upload", "Import"),
                                actionButton('reset', 'Clear'),#Clear input dataset#
                                hr(),
                                textOutput(outputId = "uploadnote"),
                                hr()),
                         column(width = 6,
                                br(),
                                # h4(strong("Data Info:")),
                                valueBoxOutput("no_cells", width = NULL),
                                valueBoxOutput("no_genes", width = NULL))))),
                fluidRow(
                  column(width = 6,
                  tabBox(
                  title = "Distinct cell types", width = NULL, id = "dc",
                  tabPanel("About", strong("Diverse cell types detected?"), textOutput(outputId = "step1_dc")),
                  tabPanel("DE", DT::DTOutput(outputId  = "dc_de_tb")%>% withSpinner(color="#FAD02C"))),
                  tabBox(
                  title = "Homogeneous cells", width = NULL, id = "homo",
                  tabPanel("About", strong("Homogeneous cells detected?"), textOutput(outputId = "step1_homogeneous")),
                  tabPanel("HVGs", DT::DTOutput(outputId  = "homo_hvgs_tb")%>% withSpinner(color="#FAD02C")),
                  tabPanel("GO",
                           helpText("Please select a database that contains information about the genome and gene annotations of a particular organism."),
                           radioButtons("go_info", "Data From",
                                        choices = list("Human" = "org.Hs.eg.db", "Mouse" = "org.Mm.eg.db"),selected = 1),
                           textOutput(outputId = "homo_go_txt"),
                           DT::DTOutput(outputId  = "homo_go_tb")%>% withSpinner(color="#FAD02C")))),
                  column(width = 6,
                         box(width=NULL, status = "warning",
                    span(textOutput(outputId = "step1_decision"), style = "font-size:20px; font-style:bold; text-align:center"),
                    plotOutput("step1_plot", height = "260px", width="550px")%>% withSpinner(color="#FAD02C")))

)),


      tabItem(tabName = "obj",
              fluidRow(
                column(3,
                      h4(strong("Generate embeddings:")),
                      "Embeddings can be generated invidivually here (in combination with the preferred trajectory method used in Step 3), or users can generate
                      their own embeddings following the workflow in the ",a(href = "https://www.rhondabacher.com/docs-escort/generate_embedding_objects_vignette.html",
                         "Generating data objects for Escort R/Shiny Vignette"), ".", 
                       numericInput("checkgenes", "The number of selected HVGs", value = 100, min=0),
                       # choose DR method
                       selectInput("checkDR", "Dimension reduction method",
                                   choices = c("MDS" = "MDS",
                                               "TSNE" = "TSNE",
                                               "UMAP" = "UMAP")),
                      br(),
                       h4(strong("Fit preliminary trajectory:")),
                       # choose Trajectory methods
                       selectInput("checkTraj", "Trajectory method", choices = c("Slingshot" = "Slingshot")),
                      br(),
                       downloadButton(outputId="downloadTraj", label = "Download .rds object")
                ),
                column(7, plotOutput("trajectory_plot")%>% withSpinner(color="#FAD02C"))


              )),

      tabItem(tabName = "step2",

              fluidRow(
                column(width=6,
                       box(
                         width=NULL, status = "primary", title = strong("Step 2: Evaluating the trajectory characteristics of embeddings"),
                         "Next, Escort identifies preferred embeddings for performing trajectory inference.
                       We will evaluate three characteristics of each embedding: ",br(),
                         " - The retention of inter-cellular relationships.",  br(),
                         " - The preservation of similarity relationships.",  br(),
                         " - Distribution of cells in the embedding space."),
                       box(
                         width=NULL, title = "Inter-cellular relationships", status = "warning",
                         tableOutput(outputId  = "step2_celltb")%>% withSpinner(color="#FAD02C")),
                       box(
                         width=NULL, title = "Preservation of similarity relationships", status = "warning",
                         tableOutput(outputId  = "step2_simitb")%>% withSpinner(color="#FAD02C")),
                       box(
                         width=NULL, title = "Cell spread ", status = "warning",
                         tableOutput(outputId = "step2_spreadtb")%>% withSpinner(color="#FAD02C"))
                       ),
                column(width=3,
                       h4(strong("Load all embeddings: (multiple allowed)")),
                       fileInput("objs", label = NULL, buttonLabel = "Upload .rds file", accept = c(".rds"), multiple = TRUE),
                       tableOutput("obj_files"))
              )
      ),


      tabItem(tabName = "step3",
              fluidRow(
                box(
                  width=6, status = "primary", title = strong("Step 3: Quantifying trajectory fitting performance"),
                  "Embeddings are also evaluated in the context of a trajectory inference method.
                A preliminary trajectory in inferred and Escort assesses the proportion of cells having an ambiguous projection
                to the trajectory. For example, trajectories in a U-shape tend to
                be inaccurate because some cells will have map to opposing pseudotimes with equal probability.")),
              fluidRow(
                box(width=6, title = "Percentage of ambiguous cells", status = "warning",
                  tableOutput(outputId = "step3_res")%>% withSpinner(color="#FAD02C")))),

      tabItem(tabName = "conclusion",
              fluidRow(
                box(
                  width=9, status = "primary", title = strong("Escort suggestions for embedding selection"),
                  "Below is a table with each embedding's rating according to their overall performance. Embeddings with a score larger than zero are
                  Recommended for trajectoy inference.",
                  hr(),
                  tableOutput(outputId = "final_res")%>% withSpinner(color="#FAD02C"))),
              downloadButton(outputId = "table", label = "Download the table"),
              hr(),
              fluidRow(
                box(
                  width=9, status = "warning",
                  plotOutput("final_plot", height = "650px", width = "900px")%>% withSpinner(color="#FAD02C")))
              )
      )
    )
  )

