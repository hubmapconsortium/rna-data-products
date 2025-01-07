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
  
  access_key_id:
    label: "AWS access key id"
    type: string
  
  secret_access_key: 
    label: "AWS secret access key"
    type: string

outputs:

  final_raw_h5ad_file:
    outputSource: secondary-analysis/final_raw_h5ad_file
    type: File
  processed_h5ad_file:
    outputSource: secondary-analysis/processed_h5ad_file
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
    outputSource: h5ad-to-arrow/output_directory
    type: Directory
  zarr_directory:
    outputSource: anndata-to-ui/output_directory
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

    out:
      - raw_h5ad_files
      - raw_h5ad_file
      - data_product_metadata
      - matrix_files
      - features_files
      - barcodes_files
    run: steps/annotate-concatenate.cwl
    label: "Annotates and concatenates h5ad data files in directory"
  
  - id: mtx-to-seurat
    scatter: [matrix_files, features_files, barcodes_files]
    scatterMethod: dotproduct
    in:
      - id: matrix_files
        source: annotate-concatenate/matrix_files
      - id: features_files
        source: annotate-concatenate/features_files
      - id: barcodes_files
        source: annotate-concatenate/barcodes_files
    
    out:
      [seurat_rds]
    run: steps/mtx-to-seurat.cwl

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
    run: steps/azimuth-annotate.cwl
    label: "Runs azimuth on the file created in the previous step"
    
  - id: add-azimuth-annotations
    in:
      - id: raw_h5ad_file
        source: annotate-concatenate/raw_h5ad_file
      - id: tissue
        source: tissue            
      - id: metadata_json
        source: azimuth-annotate/metadata_json
      - id: annotations_csv
        source: azimuth-annotate/annotations_csv
      - id: data_product_metadata
        source: annotate-concatenate/data_product_metadata

    out:
      - annotated_raw_h5ad_file
      - updated_data_product_metadata
    run: steps/add-azimuth-annotations.cwl
    
  - id: secondary-analysis
    in:
      - id: annotated_raw_h5ad_file
        source: add-azimuth-annotations/annotated_raw_h5ad_file
      - id: tissue
        source: tissue
      - id: uuids_file
        source: uuids_file
      - id: updated_data_product_metadata
        source: add-azimuth-annotations/updated_data_product_metadata
    
    out:
      - final_raw_h5ad_file
      - processed_h5ad_file
      - umap_png
      - final_data_product_metadata
    run: steps/secondary-analysis.cwl
    label: "Runs secondary anaylsis on annotated and concatenated data"

  - id: prep-to-convert
    in:
      - id: processed_h5ad_file
        source: secondary-analysis/processed_h5ad_file
    out:
      - id: secondary_analysis_dir
    run: steps/prep-to-convert.cwl
    
  - id: h5ad-to-arrow
    in:
      - id: input_directory
        source: prep-to-convert/secondary_analysis_dir
    out:
      - id: output_directory
    run: steps/h5ad-to-arrow.cwl

  - id: anndata-to-ui
    in:
      - id: input_directory
        source: prep-to-convert/secondary_analysis_dir
    out:
      - id: output_directory
    run: steps/anndata-to-ui.cwl

  - id: make-shinycell
    in:
      - id: processed_h5ad_file
        source: secondary-analysis/processed_h5ad_file
      - id: tissue
        source: tissue
      - id: metadata_file
        source: secondary-analysis/final_data_product_metadata
    out:
      - id: shinycell_dir

    run: steps/make_shinycell.cwl
    label: "Creates the shiny cell app for the data product"

  - id: upload-to-s3
    in:
      - id: final_raw_h5ad_file
        source: secondary-analysis/final_raw_h5ad_file
      - id: processed_h5ad_file
        source: secondary-analysis/processed_h5ad_file
      - id: umap_png
        source: secondary-analysis/umap_png
      - id: final_data_product_metadata
        source: secondary-analysis/final_data_product_metadata
      - id: access_key_id
        source: access_key_id
      - id: secret_access_key
        source: secret_access_key
    
    out:
      - finished_text
    
    run: steps/upload-to-s3.cwl
    label: "Uploads the pipeline outputs to s3 and ec2"

  - id: upload-zarr-to-s3
    in:
     - id: zarr_dir
       source: anndata-to-ui/output_directory
     - id: metadata_file
       source: secondary-analysis/final_data_product_metadata
     - id: access_key_id
       source: access_key_id
     - id: secret_access_key
       source: secret_access_key
    out:
     - finished_text
 
    run: steps/upload-zarr-to-s3.cwl
    label: "Uploads zarr to s3"