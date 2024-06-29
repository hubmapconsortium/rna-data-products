install.packages(
c(
  'optparse'
  ),
 Ncpus=6
)

reqPkg = c("data.table", "Matrix", "hdf5r", "reticulate", "ggplot2", 
           "gridExtra", "glue", "readr", "RColorBrewer", "R.utils", "Seurat")
newPkg = reqPkg[!(reqPkg %in% installed.packages()[,"Package"])]
if(length(newPkg)){install.packages(newPkg)}

reqPkg = c("shiny", "shinyhelper", "data.table", "Matrix", "DT", "hdf5r", 
           "ggplot2", "gridExtra", "magrittr", "ggdendro", 
           "here", "rprojroot", "tools")
newPkg = reqPkg[!(reqPkg %in% installed.packages()[,"Package"])]
if(length(newPkg)){install.packages(newPkg)}

install.packages("devtools", repos = "http://cran.rstudio.com/")

tryCatch({
    devtools::install_github("SGDDNB/ShinyCell")
},
    error = function(e) {
    message("Error installing SGDDNB/ShinyCell")
    message(e$message)
    quit("no", -1)
    }
)
