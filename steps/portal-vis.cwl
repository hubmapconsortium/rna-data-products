#!/usr/bin/env cwl-runner

class: Workflow
cwlVersion: v1.0
label: Pipeline for evaluating differential expression of genes across datasets

requirements:
  ScatterFeatureRequirement: {}

inputs:

  processed_h5ad_file:
    label: "processed_h5ad_file"
    type: File

outputs:

  arrow_directory:
    outputSource: h5ad-to-arrow/arrow_directory
    type: Directory

  zarr_store:
    outputSource: anndata-to-ui/zarr_store
    type: File
    
steps:

  - id: prep-to-convert
    in:
      - id: processed_h5ad_file
        source: secondary-analysis-pt2/processed_h5ad_file
    out:
      - id: secondary_analysis_dir
    run: portal-vis/prep-to-convert.cwl
    
  - id: h5ad-to-arrow
    in:
      - id: input_directory
        source: prep-to-convert/secondary_analysis_dir
    out:
      - id: arrow_directory
    run: portal-vis/h5ad-to-arrow.cwl

  - id: anndata-to-ui
    in:
      - id: input_directory
        source: prep-to-convert/secondary_analysis_dir
    out:
      - id: zarr_store
    run: portal-vis/anndata-to-ui.cwl