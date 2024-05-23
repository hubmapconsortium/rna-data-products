#!/usr/bin/env Rscript
library(Seurat)
library(Azimuth)
library(anndata)
library(rjson)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  cat("Usage: Rscript azimuth.R <adata_object> <ref>\n")
  quit(status = 1)
}
raw_h5ad_file <- args[1]
tissue <- args[2]

if (tissue %in% c("RK", "LK", "RL", "LL", "HT")) {
  # reference.path points to path within docker image
  if (tissue %in% c("RK", "LK")) {
    reference.path = "/opt/human_kidney"
    reference.name = "kidney"
  } else if (tissue %in% c("RL", "LL")) {
    reference.path = "/opt/human_lung"
    reference.name = "lung"
  } else if (tissue %in% "HT") { 
    reference.path = "/opt/human_heart"
    reference.name = "heart"
  } else {
    print("Tissue does not have azimuth annotations")
    # TODO: maybe put a cwl output here if i need to
    quit(save="no")
  }

  # Load reference and gather version information
  reference <- LoadReference(path = reference.path)
  reference.version <- ReferenceVersion(reference$map)
  azimuth.version <- as.character(packageVersion(pkg = "Azimuth"))
  seurat.version <- as.character(packageVersion(pkg = "Seurat"))
  max.dims <- as.double(length(slot(reference$map, "reductions")$refDR))
  meta.data <- names(slot(reference$map, "meta.data"))
  annotation.columns <- c()
  if (reference.name %in% c("lung")) {
    annotation.columns <- c("ann_level_1", "ann_level_2", "ann_level_3", "ann_level_4", "ann_level_5", "ann_finest_level")
  } else {
    for (i in grep("[.]l[1-3]", meta.data)) {
      annotation.columns <- c(annotation.columns, meta.data[i])
    }
  }

  # Load the query object for mapping
  # Change the file path based on where the query file is located on your system.
  query <- LoadFileInput(path = raw_h5ad_file)
  cat("Query object dimensions: ", dim(query), "\n")
  raw_h5ad <- read_h5ad(filename = raw_h5ad_file)
  cat("Raw h5ad object dimensions: ", dim(raw_h5ad), "\n")
  query <- ConvertGeneNames(
    object = query,
    reference.names = rownames(x = reference$map),
    homolog.table = 'https://seurat.nygenome.org/azimuth/references/homologs.rds'
  )
  cat("Query object dimensions: ", dim(query), "\n")
  # Calculate nCount_RNA and nFeature_RNA if the query does not
  # contain them already
  if (!all(c("nCount_RNA", "nFeature_RNA") %in% c(colnames(x = query[[]])))) {
      calcn <- as.data.frame(x = Seurat:::CalcN(object = query))
      colnames(x = calcn) <- paste(
        colnames(x = calcn),
        "RNA",
        sep = '_'
      )
      query <- AddMetaData(
        object = query,
        metadata = calcn
      )
      rm(calcn)
  }
  cat("Query object dimensions: ", dim(query), "\n")
  # Calculate percent mitochondrial genes if the query contains genes
  # matching the regular expression "^MT-"
  if (any(grepl(pattern = '^MT-', x = rownames(x = query)))) {
    query <- PercentageFeatureSet(
      object = query,
      pattern = '^MT-',
      col.name = 'percent.mt',
      assay = "RNA"
    )
  }
  cat("Query object dimensions: ", dim(query), "\n")
  assay <- DefaultAssay(query)
  anchors <- FindTransferAnchors(
    reference = reference$map,
    query = query,
    k.filter = NA,
    reference.neighbors = "refdr.annoy.neighbors",
    reference.assay = "refAssay",
    query.assay = assay, 
    reference.reduction = "refDR",
    normalization.method = "SCT",
    features = rownames(Loadings(reference$map[["refDR"]])),
    dims = 1:max.dims,
    n.trees = 20,
    mapping.score.k = 100
  )
  cat("Query object dimensions: ", dim(query), "\n")

  # Fix for NA values 
  cell_embeddings <- anchors@object.list[[1]][["pcaproject"]]@cell.embeddings
  cell_embeddings[rowSums(is.na(cell_embeddings)) == ncol(cell_embeddings), ] <- 0
  anchors@object.list[[1]][["pcaproject"]]@cell.embeddings <- cell_embeddings
  cat("Query object dimensions: ", dim(query), "\n")
  # Transfer cell type labels and impute protein expression
  #
  # Transferred labels are in metadata columns named "predicted.*"
  # The maximum prediction score is in a metadata column named "predicted.*.score"
  # The prediction scores for each class are in an assay named "prediction.score.*"
  # The imputed assay is named "impADT" if computed

  refdata <- lapply(X = annotation.columns, function(x) {
    reference$map[[x, drop = TRUE]]
  })
  names(x = refdata) <- annotation.columns

  query <- TransferData(
    reference = reference$map,
    query = query,
    dims = 1:max.dims,
    anchorset = anchors,
    refdata = refdata,
    n.trees = 20,
    store.weights = TRUE
  )
  cat("Query object dimensions: ", dim(query), "\n")
  # Calculate the embeddings of the query data on the reference SPCA
  query <- IntegrateEmbeddings(
    anchorset = anchors,
    reference = reference$map,
    query = query,
    reductions = "pcaproject",
    reuse.weights.matrix = TRUE
  )
  cat("Query object dimensions: ", dim(query), "\n")
  # Calculate the query neighbors in the reference
  # with respect to the integrated embeddings
  query[["query_ref.nn"]] <- FindNeighbors(
    object = Embeddings(reference$map[["refDR"]]),
    query = Embeddings(query[["integrated_dr"]]),
    return.neighbor = TRUE,
    l2.norm = TRUE
  )
  cat("Query object dimensions: ", dim(query), "\n")
  # The reference used in the app is downsampled compared to the reference on which
  # the UMAP model was computed. This step, using the helper function NNTransform,
  # corrects the Neighbors to account for the downsampling.
  query <- Azimuth:::NNTransform(
    object = query,
    meta.data = reference$map[[]]
  )
  cat("Query object dimensions: ", dim(query), "\n")
  # Project the query to the reference UMAP.
  query[["proj.umap"]] <- RunUMAP(
    object = query[["query_ref.nn"]],
    reduction.model = reference$map[["refUMAP"]],
    reduction.key = 'UMAP_'
  )
  cat("Query object dimensions: ", dim(query), "\n")
  umap_embeddings <- query[["proj.umap"]]@cell.embeddings
  
  # Calculate mapping score and add to metadata
  query <- AddMetaData(
    object = query,
    metadata = MappingScore(anchors = anchors, ndim = max.dims),
    col.name = "mapping.score"
  )
  cat("Query object dimensions: ", dim(query), "\n")
  # build and save df containing annotations and scores
  # need to gather column names to save based on the names of things which are in the reference
  # we know which columns exist based on which reference is used, thus include the columns and
  # annotation levels as metadata
  predicted.cols <- c()
  for (col in annotation.columns) {
    predicted.cols <- c(predicted.cols, paste0("predicted.", col), paste0("predicted.", col, ".score"))
  }

  # build list of matrices to append to anndata object
  ls <- list()
  for (col in predicted.cols) {
    ls <-  c(ls, list(matrix(query[[col]][[1]])))
  }
  names(ls) <- predicted.cols
  ls <- c(ls, list("mapping.score" = matrix(query$mapping.score)))

  # build and save data as a CSV
  df <- data.frame(ls)
  df <- cbind(cells = rownames(df), df)
  df["barcodes"] <- names(query$orig.ident)
  df["V1"] <- matrix(query[["proj.umap"]]@cell.embeddings, ncol=2)[,1]
  df["V2"] <- matrix(query[["proj.umap"]]@cell.embeddings, ncol=2)[,2]

  # load secondary analysis matrix as anndata object
  raw_h5ad <- read_h5ad(filename = raw_h5ad_file)

  # add reference-guided UMAP to anndata object
  umap.new <- umap_embeddings
  raw_h5ad$obsm$X_umap_proj <- umap.new

  # save modified secondary_analysis.h5ad matrix to a new annotated equivalent
  write_h5ad(raw_h5ad, secondary.analysis.path) 
  write.csv(df, file=annotations.csv.path, row.names=FALSE)

  version.metadata <- list(
    "is_annotated" = TRUE,
    "seurat" = list("version" = seurat.version),
    "azimuth" = list("version" = azimuth.version),
    "azimuth_reference" = list("version" = reference.version, "name" = reference.name),
    "annotation_names" = predicted.cols
  )
  version.metadata.json = toJSON(version.metadata)
  f <- file("version_metadata.json")
  write(version.metadata.json, f)
  close(f)

} else {
  # no-op, but still return the unmodified secondary_analysis.h5ad and metadata indicating no annotation occurred
  write.csv(data.frame(), file=annotations.csv.path, row.names=FALSE)
  ad <- read_h5ad(raw_h5ad)
  write_h5ad(ad, raw_h5ad)
  version.metadata <- list("is_annotated" = FALSE)
  version.metadata.json = toJSON(version.metadata)
  f <- file(version.metadata.path)
  write(version.metadata.json, f)
  close(f)
  # Create dummy annotations file if no annotation performed. Will handle this case in write_metadata.py
}
