install.packages(
c(
  'optparse'
  ),
 Ncpus=6
)

install.packages("BiocManager", repos = "https://cloud.r-project.org")
# bio_pkgs <- c("Biostrings", "BSgenome", "CNEr", "DirichletMultinomial", "GenomeInfoDb", "GenomicRanges", "IRanges", "rtracklayer", "S4Vectors", "XVector", "glmGamPoi", "TFBSTools", "S4Arrays",
#                 "SparseArray", "DelayedArray", "")
# tryCatch({
#     BiocManager::install("DirichletMultinomial")
# })

tryCatch({
    install.packages("magick", version = "2.7.3")
},
    error = function(e) {
    message("Error installing magick")
    message(e$message)
    quit("no", -1)
  }
)

tryCatch({
    devtools::install_github("satijalab/seurat", "seurat5")
},
    error = function(e) {
    message("Error installing satijalab/seurat")
    message(e$message)
    quit("no", -1)
    }
)

tryCatch({
    devtools::install_github("satijalab/seurat-data", "seurat5")
},
    error = function(e) {
    message("Error installing satijalab/seurat-data")
    message(e$message)
    quit("no", -1)
    }
)

tryCatch({
    devtools::install_github("satijalab/azimuth", "master")

},
    error = function(e) {
    message("Error installing satijalab/azimuth")
    message(e$message)
    quit("no", -1)
    }
)
reqPkg = c("data.table", "Matrix", "hdf5r", "reticulate", "ggplot2", 
           "gridExtra", "glue", "readr", "RColorBrewer", "R.utils", "Seurat")
newPkg = reqPkg[!(reqPkg %in% installed.packages()[,"Package"])]
if(length(newPkg)){install.packages(newPkg)}

tryCatch({
    reticulate::py_install("anndata")
},
    error = function(e) {
    message("Error installing anndata")
    message(e$message)
    quit("no", -1)
    }
)

reqPkg = c("shiny", "shinyhelper", "data.table", "Matrix", "DT", "hdf5r", 
           "reticulate", "ggplot2", "gridExtra", "magrittr", "ggdendro")
newPkg = reqPkg[!(reqPkg %in% installed.packages()[,"Package"])]
if(length(newPkg)){install.packages(newPkg)}

tryCatch({
    devtools::install_github("SGDDNB/ShinyCell")
},
    error = function(e) {
    message("Error installing SGDDNB/ShinyCell")
    message(e$message)
    quit("no", -1)
    }
)
