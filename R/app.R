# Run the application
#' Escort
#'
#' Runs the Escort Shiny app.
#'
#' @return This starts up the shiny app interface for Escort.
#' @importFrom shinyjs show alert
#' @import mclust
#' @importFrom stats predict
#' @importFrom scales alpha
#' @importFrom dplyr first setequal intersect union 
#' @import magrittr
#' @importFrom shinycssloaders withSpinner
#' @import DT
#' @import shinydashboard
#' @import parallelDist
#' @import org.Hs.eg.db
#' @import org.Mm.eg.db
#' @export
#'

shinyEscort <- function() {
  options(
    shiny.maxRequestSize = 10000*1024^5,
    shiny.launch.browser = FALSE,
    shiny.port = 19901,
    shiny.host = "0.0.0.0",
    test.mode = getOption("shiny.testmode", FALSE),
    )
  library(Escort)
  shinyApp(ui, server)
}

