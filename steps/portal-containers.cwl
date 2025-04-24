#!/usr/bin/env cwl-runner

class: Workflow
cwlVersion: v1.0
label: Pipeline for evaluating differential expression of genes across datasets

requirements:
  ScatterFeatureRequirement: {}

inputs:

  processed_h5ad_file:
    type: File

outputs:

  arrow_directory:
    outputSource: h5ad-to-arrow/output_directory
    type: Directory
  zarr_directory:
    outputSource: anndata-to-ui/output_directory
    type: Directory
    
steps:

  - id: prep-to-convert
    in:
      - id: processed_h5ad_file
        source: secondary-analysis-pt2/processed_h5ad_file
    out:
      - id: secondary_analysis_dir
    run: steps/portal-containers/prep-to-convert.cwl
    
  - id: h5ad-to-arrow
    in:
      - id: input_directory
        source: prep-to-convert/secondary_analysis_dir
    out:
      - id: output_directory
    run: steps/portal-containers/h5ad-to-arrow.cwl

  - id: anndata-to-ui
    in:
      - id: input_directory
        source: prep-to-convert/secondary_analysis_dir
    out:
      - id: output_directory
    run: steps/portal-containers/anndata-to-ui.cwl