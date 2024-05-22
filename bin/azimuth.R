#!/usr/bin/env Rscript
library(Seurat)
library(Azimuth)
library(SeuratData)
library(patchwork)
library(SeuratDisk)
library(magick)
library(R.utils)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  cat("Usage: Rscript azimuth.R <adata_object> <ref>\n")
  quit(status = 1)
}
raw_h5ad_file <- args[1]
tissue <- args[2]

# wd <- getwd()
# copyDirectory(wd, "azimuth-output")
# file.copy(raw_h5ad_file, "azimuth_output")
# raw_h5ad_path <- sprintf("azimuth_output/%", raw_h5ad_file)

if (tissue %in% c("RK", "LK", "RL", "LL", "HT")) {
  if (tissue %in% c("RK", "LK")) {
    ref = "kidneyref"
  } else if (tissue %in% c("RL", "LL")) {
    ref = "lungref"
  } else if (tissue %in% "HT") { 
    ref = "heartref"
  } else {
    print("Tissue does not have azimuth annotations")
    # TODO: maybe put a cwl output here if i need to
    quit(save="no")
  }
}
print("Running Azimuth")
azimuth_object <-RunAzimuth(query=raw_h5ad_file, reference=ref)
write.csv(azimuth_object@meta.data, "azimuth_results.csv")
