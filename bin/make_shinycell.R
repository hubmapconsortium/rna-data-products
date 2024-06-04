library(ShinyCell)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  cat("Usage: Rscript make_shinycell.R <adata_object> <ref>\n")
  quit(status = 1)
}
inpFile <- args[1]
tissue <- args[2]


scConf = createConfig(inpFile)
#dir = sprintf("shinyApps/% s/", tissue)
dir = "shinyApps/"
title = sprintf("Shiny Cell h5ad % s", tissue)
makeShinyApp(inpFile, scConf, shiny.dir = dir, shiny.title = title) 
