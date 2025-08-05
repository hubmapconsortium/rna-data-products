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
  
  organism:
    type: string?
    default: "human"

outputs:

  final_raw_h5mu_file:
    outputSource: secondary-analysis/final_raw_h5mu_file
    type: File
  processed_h5mu_file:
    outputSource: secondary-analysis/processed_h5mu_file
    type: File
  umap_png:
    outputSource: secondary-analysis/umap_png
    type: File
  final_data_product_metadata:
    outputSource: secondary-analysis/final_data_product_metadata
    type: File
  shinycell_dir:
    outputSource: make-shinycell/shinycell_dir
    type: Directory
  arrow_directory:
    outputSource: portal-vis/arrow_directory
    type: Directory
  zarr_store:
    outputSource: portal-vis/zarr_store
    type: Directory
    
steps:

  - id: annotate-concatenate
    in:
      - id: enable_manhole
        source: enable_manhole
      - id: data_directory
        source: data_directory
      - id: uuids_file
        source: uuids_file
      - id: tissue
        source: tissue
      - id: organism
        source: organism

    out:
      - annotated_raw_h5ad_file
      - updated_data_product_metadata
    run: steps/annotate-concatenate.cwl
    label: "Annotates and concatenates h5ad data files in directory"

  - id: secondary-analysis
    in:
      - id: annotated_raw_h5ad_file
        source: annotate-concatenate/annotated_raw_h5ad_file
      - id: tissue
        source: tissue
      - id: uuids_file
        source: uuids_file
      - id: updated_data_product_metadata
        source: annotate-concatenate/updated_data_product_metadata
    
    out:
      - final_raw_h5mu_file
      - processed_h5ad_file
      - processed_h5mu_file
      - umap_png
      - final_data_product_metadata
    run: steps/secondary-analysis.cwl
    label: "Runs secondary anaylsis on annotated and concatenated data"
  
  - id: portal-vis
    in:
      - id: processed_h5ad_file
        source: secondary-analysis/processed_h5ad_file
    out:
      - arrow_directory
      - zarr_store
    run: steps/portal-vis.cwl

  - id: make-shinycell
    in:
      - id: processed_h5ad_file
        source: secondary-analysis/processed_h5ad_file
      - id: tissue
        source: tissue
      - id: metadata_file
        source: secondary-analysis/final_data_product_metadata
    out:
      - shinycell_dir

    run: steps/make_shinycell.cwl
    label: "Creates the shiny cell app for the data product"

  - id: upload-to-vm
    in:
      - id: final_raw_h5mu_file
        source: secondary-analysis/final_raw_h5mu_file
      - id: processed_h5mu_file
        source: secondary-analysis/processed_h5mu_file
      - id: umap_png
        source: secondary-analysis/umap_png
      - id: final_data_product_metadata
        source: secondary-analysis/final_data_product_metadata
      - id: zarr_store
        source: portal-vis/zarr_store
    
    out:
      - finished_text
    
    run: steps/upload-to-vm.cwl
    label: "Uploads the pipeline outputs to vm"