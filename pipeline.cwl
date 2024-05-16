#!/usr/bin/env cwl-runner

class: Workflow
cwlVersion: v1.0
label: Pipeline for evaluating differential expression of genes across datasets

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
  
  access_key_id:
    label: "AWS access key id"
    type: string
  
  secret_access_key: 
    label: "AWS secret access key"
    type: string

outputs:

  raw_h5ad_file:
    outputSource: secondary-analysis/raw_h5ad_file
    type: File
  processed_h5ad_file:
    outputSource: secondary-analysis/processed_h5ad_file
    type: File
  umap_png:
    outputSource: secondary-analysis/umap_png
    type: File

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

    out:
      - raw_h5ad_file
    run: steps/annotate-concatenate.cwl
    label: "Annotates and concatenates h5ad data files in directory"
  
  - id: azimuth-annotate
    in: 
      - id: raw_h5ad_file
        source: annotate-concatenate/raw_h5ad_file
      - id: tissue:
        source: tissue
    
    out:
      - annotations_csv
    run: steps/azimuth-annotate
    label: "Runs azimuth on the file created in the previous step"

  - id: secondary-analysis
    in:
      - id: raw_h5ad_file
        source: annotate-concatenate/raw_h5ad_file
      - id: annotated_csv
        source: azimuth-annotate/annotations_csv
      - id: tissue
        source: tissue
    
    out:
      - raw_h5ad_file
      - processed_h5ad_file
      - umap_png
    run: steps/secondary-analysis.cwl
    label: "Runs secondary anaylsis on annotated and concatenated data"

  - id: upload-to-s3
    in:
      - id: raw_h5ad_file
        source: secondary-analysis/annotated_raw_h5ad_file
      - id: processed_h5ad_file
        source: secondary-analysis/processed_h5ad_file
      - id: umap_png
        source: secondary-analysis/umap_png
      - id: tissue
        source: tissue
      - id: access_key_id
        source: access_key_id
      - id: secret_access_key
        source: secret_access_key
    
    out:
      - finished_text
    
    run: steps/upload-to-s3.cwl
    label: "Uploads the pipeline outputs to s3"