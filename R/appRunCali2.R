#' @export
cali_app2 <- function() {
  appDir <- system.file("cali2",
                        package = "spatioTemp")

  shiny::runApp(
    appDir,
    launch.browser = TRUE,
    display.mode = "normal"
  )
}
