library(ShinyCell)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  cat("Usage: Rscript make_shinycell.R <adata_object> <ref>\n")
  quit(status = 1)
}
inpFile <- args[1]
tissue <- args[2]


scConf = createConfig(inpFile)
mainDir = "shinyApps"
subDir = tissue
shinyDir <- file.path(mainDir, subDir)
if (!dir.exists(mainDir)) {dir.create(mainDir)}
if (!dir.exists(shinyDir)) {dir.create(shinyDir)}  

title = sprintf("Shiny Cell h5ad % s", tissue)
makeShinyApp(inpFile, scConf, shiny.dir = shinyDir, shiny.title = title) 
