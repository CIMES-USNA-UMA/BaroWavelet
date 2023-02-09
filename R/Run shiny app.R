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
#' @examples
#' RunBaroWaveletApp()
#' 
RunBaroWaveletApp <- function(){
  shiny::runGitHub(repo = "CIMES-USNA-UMA/BaroWaveletApp",
                launch.browser = TRUE)
}
