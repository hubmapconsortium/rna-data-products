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
adata_object <- args[1]
tissue <- args[2]

getwd()

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
azimuth_object <-RunAzimuth(query=adata_object, reference=ref)
write.csv(azimuth_object@meta.data, "azimuth_results.csv")
