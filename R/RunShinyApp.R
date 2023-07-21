#' Run BaroWavelet's ShinyApp
#'
#' Launches BaroWaveletApp, BaroWavelet's ShinyApp
#'
#'
#' @author Alvaro Chao-Ecija
#' 
#'
#' @export
#'
#'
#' @examples
#' \dontrun{
#' RunBaroWaveletApp()
#' }
#'
RunBaroWaveletApp <- function() {
  if (!requiredNamespace("shiny"))
    stop(
      "Package 'shiny' is not installed. Package 'shiny' must be installed to 
      access the shiny aplication."
    )
  shiny::runGitHub(repo = "CIMES-USNA-UMA/BaroWaveletApp",
                   launch.browser = TRUE)
}
