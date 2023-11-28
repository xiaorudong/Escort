# Run the application
#' Escort
#'
#' Runs the Escort Shiny app.
#'
#' @return This starts up the shiny app interface for Escort.
#' @import Seurat
#' @import shiny
#' @import shinyWidgets
#' @import mclust
#' @import slingshot
#' @importFrom scales alpha
#' @import RColorBrewer
#' @import dplyr
#' @importFrom shinycssloaders withSpinner
#' @import DT
#' @import shinydashboard
#' @import shinyjs
#' @import parallelDist
#' @import org.Hs.eg.db
#' @import org.Mm.eg.db
#' @export
#'

shinyEscort <- function() {
  options(shiny.maxRequestSize = 10000*1024^5)
  shinyApp(ui, server)
}


