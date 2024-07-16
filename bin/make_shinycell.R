library(ShinyCell)
library(rjson)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  cat("Usage: Rscript make_shinycell.R <adata_object> <ref>\n")
  quit(status = 1)
}
inpFile <- args[1]
tissue <- args[2]
metadata_file <- args[3]

tryCatch({
  scConf = createConfig(inpFile)
},
error = function(e) {
  reticulate::py_last_error()
  message("Error creating config")
  message(e$message)
  quit("no", -1)
}
)
mainDir = "shinyApps"
subDir = tissue

json_data <- fromJSON(file=metadata_file)
uuid = json_data['Data Product UUID']

tissueDir <- file.path(mainDir, subDir)
shinyDir <- file.path(mainDir, subDir, uuid)
if (!dir.exists(mainDir)) {dir.create(mainDir)}
if (!dir.exists(tissueDir)) {dir.create(tissueDir)}
if (!dir.exists(shinyDir)) {dir.create(shinyDir)}  

title = sprintf("Shiny Cell h5ad % s", tissue)
makeShinyApp(inpFile, scConf, shiny.dir = shinyDir, shiny.title = title) 
