
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
    sidebarMenu(id="tabs",
                menuItem("About", tabName = "home", icon = icon("house-user")),
                menuItem("Step 1", tabName = "step1"),
                menuItem("Downstrem Analysis", tabName = "obj"),
                menuItem("Step 2", tabName = "step2"),
                menuItem("Step 3", tabName = "step3"),
                menuItem("Conclusion", tabName = "conclusion")
                #
                # menuItem("Start Analysis", tabName = "app", icon = icon("circle-play"),
                #          startExpanded = TRUE,
                #          # menuSubItem("Load Data", tabName = "load"),
                #          menuSubItem("Step 1", tabName = "step1"),
                #          menuSubItem("Downstrem Analysis", tabName = "obj"),
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
              h3(strong("Welcome to Escort")),
              fluidRow(
                column(6,
                       "Escort is a framework that uses data to evaluate
                       how well different methods can infer trajectories
                       from single-cell data. It employs a three-step approach to guide users
                       through the analysis by providing evaluations of embeddings
                       representing a series of analysis choices including feature selection,
                       dimension reduction, and trajectory inference methods.",
                       br(),
                       br(),
                       "The first step is to check if the data have a trajectory signal or not.
                       Sometimes, the data are not suitable for trajectory analysis,
                       for example, when the cells come from biologically distinct clusters or
                       sets of cells with insufficient heterogeneity. In these cases,
                       Escort will tell the user that trajectory inference is not a good idea.",
                       br(),
                       br(),
                       "The second step is to compare different low-dimensional embeddings and
                       see how well they capture the relationships between cells
                       in the original high-dimensional space.
                       Escort also looks at how dense the cells are in the low-dimensional space,
                       and prefers embeddings that have a tighter distribution for trajectory analysis.",
                       br(),
                       br(),
                       "The third step is to evaluate how well a specific trajectory inference method
                       works with a given embedding. This depends on the graph structure and other parameters
                       that the method uses. Escort also measures how clear the projection of each cell
                       is along the trajectory, and avoids methods that produce ambiguous results.
                       Finally, Escort gives an overall score for each analysis option based on the three steps,
                       and helps the user to find the best way to infer a trajectory from their data. Go to",
                       a(href = "https://github.com/NabiilahArdini/Shiny-Box",
                         "GitHub Page"),
                       "to find more details on the source code.", br()),
                column(6, imageOutput("home_img"))
                )),

      tabItem(tabName = "step1",
              fluidRow(
                column(width=6,
                       h4(strong("Upload scRNA-seq datasets:")),
                       "Note: please upload csv files for raw data and normalized data.")),
              fluidRow(column(width = 3,
                              # read raw data after QC
                              fileInput("rawfile", label = "Raw Data", buttonLabel = "Upload", accept = c(".csv")),
                              # read norm data after QC
                              fileInput("normfile", label = "Norm Data", buttonLabel = "Upload", accept = c(".csv")),
                              actionButton("upload", "Import"),
                              actionButton('reset', 'Clear'),#Clear input dataset#
                              hr(),
                              textOutput(outputId = "uploadnote"),
                              hr()),
                       column(width = 3,
                              # h4(strong("Data Info:")),
                              valueBoxOutput("no_cells", width = NULL),
                              valueBoxOutput("no_genes", width = NULL))),

              fluidRow(
                column(width = 6,
                       box(width=NULL, status = "primary", title = strong("Step 1: Detecting Trajectory Existence"),
                           "In the first step of our analysis, we will identify two scenarios
                where trajectory fitting is not feasible: ",br(),
                           " - Cells represent diverse cell types",  br(),
                           " - Cells are homogeneous"),

                       tabBox(
                         title = "Diverse cell types", width = NULL, id = "dc",
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
                                  DT::DTOutput(outputId  = "homo_go_tb")%>% withSpinner(color="#FAD02C")))

                       ),

                column(width = 6,
                       box(width=NULL, status = "warning",
                           span(textOutput(outputId = "step1_decision"), style = "font-size:20px; font-style:bold; text-align:center"),
                           plotOutput("step1_plot", height = "260px", width="550px")%>% withSpinner(color="#FAD02C"))
                ))),


      tabItem(tabName = "obj",
              fluidRow(
                column(3,
                      h4(strong("Generate embeddings:")),
                       numericInput("checkgenes", "The number of selected HVGs", value = 100, min=0),
                       # choose DR method
                       selectInput("checkDR", "Dimension Reduction Method",
                                   choices = c("MDS" = "MDS",
                                               "TSNE" = "TSNE",
                                               "UMAP" = "UMAP")),
                      br(),
                       h4(strong("Fit Trajectory:")),
                       # choose Trajectory methods
                       selectInput("checkTraj", "Trajectory Methods", choices = c("Slingshot" = "Slingshot")),
                      br(),
                       downloadButton(outputId="downloadTraj", label = "Download Trajectory")
                ),
                column(7, plotOutput("trajectory_plot")%>% withSpinner(color="#FAD02C"))


              )),

      tabItem(tabName = "step2",

              fluidRow(
                column(width=6,
                       box(
                         width=NULL, status = "primary", title = strong("Step 2: Evaluating the Characteristics of Embeddings"),
                         "The second step is designed to identify preferred embeddings for performing trajectory inference.
                Since all methods employ some form of dimension reduction,
                the first evaluation is the efficacy of low-dimensional embeddings in maintaining
                the inter-cellular relationships found in the high-dimensional data.
                The accuracy of trajectory prediction is heavily dependent on the extent to
                which these relationships are preserved in the embedding.
                We will evaluate three characteristics of embeddings: ",br(),
                         " - The retention of inter-cellular relationships that are present in the high-dimensional data",  br(),
                         " - The preservation of similarity relationships in the low-dimensional embedding",  br(),
                         " - Cell spread in embeddings"),
                       box(
                         width=NULL, title = "Inter-cellular Relationships", status = "warning",
                         tableOutput(outputId  = "step2_celltb")%>% withSpinner(color="#FAD02C")),
                       box(
                         width=NULL, title = "Preservation of Similarity Relationships", status = "warning",
                         tableOutput(outputId  = "step2_simitb")%>% withSpinner(color="#FAD02C")),
                       box(
                         width=NULL, title = "Cell Spread ", status = "warning",
                         tableOutput(outputId = "step2_spreadtb")%>% withSpinner(color="#FAD02C"))
                       ),
                column(width=3,
                       h4(strong("Load Embeddings:")),
                       fileInput("objs", label = NULL, buttonLabel = "Upload rds file", accept = c(".rds"), multiple = TRUE),
                       tableOutput("obj_files"))
              )
      ),


      tabItem(tabName = "step3",
              fluidRow(
                box(
                  width=6, status = "primary", title = strong("Step 3: Quantifying trajectory fitting performance"),
                  "Now that the embeddings have been evaluated independently,
                I next evaluate them in the context of a trajectory inference method.
                When using a particular method to fit a trajectory, Escort assessed
                how many cells are positioned along the trajectory such that their projection
                is ambiguous. For example, trajectories in a U-shape tend to
                be inaccurate because some cells can be mapped with similar probabilities to
                either the starting or ending points of the trajectory.")),
              fluidRow(
                box(width=6, title = "Percentage of Ambiguous Cells", status = "warning",
                  tableOutput(outputId = "step3_res")%>% withSpinner(color="#FAD02C")))),

      tabItem(tabName = "conclusion",
              fluidRow(
                box(
                  width=9, status = "primary", title = strong("Suggestions for Embedding Selection"),
                  "Given the above analysis, we provide a table showing the recommendations for
                embedding choices based on whether embedding can preserve the original
                data structure and whether the trajectory adequately reflects
                the differentiation of the cells. We rate each embedding
                according to their overall performance and provide suggestions
                for selecting them.",
                  hr(),
                  tableOutput(outputId = "final_res")%>% withSpinner(color="#FAD02C"))),
              fluidRow(
                box(
                  width=9, status = "warning",
                  plotOutput("final_plot", height = "650px", width = "900px")%>% withSpinner(color="#FAD02C")))
              )
      )
    )
  )

