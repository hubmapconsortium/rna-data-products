#!/usr/bin/env cwl-runner

class: Workflow
cwlVersion: v1.0
label: Pipeline for evaluating differential expression of genes across datasets

requirements:
  ScatterFeatureRequirement: {}

inputs:

  enable_manhole:
    label: "Whether to enable remote debugging via 'manhole'"
    type: boolean?

  data_directory:
    label: "Path to directory containing processed RNA datasets"
    type: Directory

  uuids_file:
    label: "Path to a file containing a list of uuids for the dataset to be indexed"
    type: File

  tissue:
    label: "String description of tissue type"
    type: string?

outputs:

  annotated_raw_h5ad_file:
    outputSource: add-azimuth-annotations/annotated_raw_h5ad_file
    type: File
  updated_data_product_metadata:
    outputSource: add-azimuth-annotations/updated_data_product_metadata
    type: File
    
steps:

  - id: concatenate
    in:
      - id: enable_manhole
        source: enable_manhole
      - id: data_directory
        source: data_directory
      - id: uuids_file
        source: uuids_file
      - id: tissue
        source: tissue

    out:
      - raw_h5ad_files
      - raw_h5ad_file
      - data_product_metadata
      - matrix_files
      - features_files
      - barcodes_files
    run: annotate-concatenate/concatenate.cwl
    label: "Concatenates h5ad data files in directory"
  
  - id: mtx-to-seurat
    scatter: [matrix_files, features_files, barcodes_files]
    scatterMethod: dotproduct
    in:
      - id: matrix_files
        source: concatenate/matrix_files
      - id: features_files
        source: concatenate/features_files
      - id: barcodes_files
        source: concatenate/barcodes_files
    
    out:
      [seurat_rds]
    run: annotate-concatenate/mtx-to-seurat.cwl

  - id: azimuth-annotate
    scatter: [seurat_rds]
    scatterMethod: dotproduct
    in: 
      - id: seurat_rds
        source: mtx-to-seurat/seurat_rds
      - id: tissue
        source: tissue
    
    out:
      - annotations_csv
      - metadata_json
    run: annotate-concatenate/azimuth-annotate.cwl
    label: "Runs azimuth on the file created in the previous step"
    
  - id: add-azimuth-annotations
    in:
      - id: raw_h5ad_file
        source: concatenate/raw_h5ad_file
      - id: tissue
        source: tissue            
      - id: metadata_json
        source: azimuth-annotate/metadata_json
      - id: annotations_csv
        source: azimuth-annotate/annotations_csv
      - id: data_product_metadata
        source: concatenate/data_product_metadata

    out:
      - annotated_raw_h5ad_file
      - updated_data_product_metadata
    run: annotate-concatenate/add-azimuth-annotations.cwl