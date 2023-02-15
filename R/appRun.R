#' @export
valencia_app <- function() {
  appDir <- system.file("valencia",
                        package = "spatioTemp")

  shiny::runApp(
    appDir,
    launch.browser = TRUE,
    display.mode = "normal"
  )
}



