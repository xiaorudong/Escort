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

shinyEscort <- function(port = 19901, is_local_development = TRUE) {
  # default port is set to 19901
  # if is_local_development is set to false, then a browser window will
  # not open up. We should use FALSE in production server
  options_args <- list(
    shiny.maxRequestSize = 10000 * 1024^5,
    shiny.launch.browser = is_local_development,
    shiny.port = port,
    test.mode = getOption("shiny.testmode", FALSE)
  )
  if (is_local_development == FALSE) {
    options_args$shiny.host <- "0.0.0.0"
  }
  do.call(options, options_args)
  library(Escort)
  shinyApp(ui, server)
}
