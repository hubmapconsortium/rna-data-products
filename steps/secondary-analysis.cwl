#!/usr/bin/env cwl-runner

class: Workflow
cwlVersion: v1.0
label: Pipeline for evaluating differential expression of genes across datasets

requirements:
  ScatterFeatureRequirement: {}

inputs:

  annotated_raw_h5ad_file:
    type: File

  updated_data_product_metadata:
    label: "Path to directory containing processed RNA datasets"
    type: File

  uuids_file:
    label: "Path to a file containing a list of uuids for the dataset to be indexed"
    type: File

  tissue:
    label: "String description of tissue type"
    type: string?    

outputs:

  annotated_raw_h5ad_file:
    outputSource: add-azimuth-annotations/final_raw_h5ad_file
    type: File
  updated_data_product_metadata:
    outputSource: add-azimuth-annotations/final_data_product_metadata
    type: File
    
steps:

  - id: secondary-analysis-pt1
    in:
      - id: annotated_raw_h5ad_file
        source: annotated_raw_h5ad_file
      - id: tissue
        source: tissue
      - id: uuids_file
        source: uuids_file
      - id: updated_data_product_metadata
        source: updated_data_product_metadata
    
    out:
      - final_raw_h5ad_file
      - partially_processed_h5ad_file
      - updated_data_product_metadata
    run: steps/seconadry-analysis/secondary-analysis-pt1.cwl
    label: "Runs secondary anaylsis on annotated and concatenated data"
  
  - id: sketching
    in:
      - id: partially_processed_h5ad_file
        source: secondary-analysis-pt1/partially_processed_h5ad_file
    
    out:
      - sketched_h5ad_file
    run: steps/secondary-analysis/sketching.cwl
  
  - id: secondary-analysis-pt2
    in:
      - id: sketched_h5ad_file
        source: sketching/sketched_h5ad_file
      - id: tissue
        source: tissue
      - id: updated_data_product_metadata
        source: updated_data_product_metadata
    
    out:
      - processed_h5ad_file
      - umap_png
      - final_data_product_metadata
    run: steps/secondary-analysis/secondary-analysis-pt2.cwl
    label: "Runs secondary anaylsis on annotated and concatenated data"
