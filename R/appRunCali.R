#' @export
cali_app <- function() {
  appDir <- system.file("cali",
                        package = "spatioTemp")

  shiny::runApp(
    appDir,
    launch.browser = TRUE,
    display.mode = "normal"
  )
}
