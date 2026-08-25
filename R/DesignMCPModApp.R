#' Start externally hosted DesignMCPMod Shiny App
#' 
#' 
#' This function starts the externally hosted DesignMCPMod Shiny App in a
#' browser window. The app was developed by Sophie Sun \[aut, cre\], Danyi Xiong
#' \[aut\], Bjoern Bornkamp \[ctb\], Frank Bretz \[ctb\], Ardalan Mirshani \[ctb\].
#' This app performs power and sample size calculations for a multiple contrast
#' test for normal, binary and negative binomial outcomes. The app uses the
#' DoseFinding package as calculation backend and the R code underlying the
#' calculations in the app can be extracted from the app.
#' 
#' @export
DesignMCPModApp <- function(){
  browseURL("https://huisophiesunrshiny.shinyapps.io/designmcpmod/")
}
