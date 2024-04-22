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

outputs:

  raw_h5ad_file:
    outputSource: annotate-concatenate/raw_h5ad_file
    type: File
  processed_h5ad_file:
    outputSource: annotate-concatenate/processed_h5ad_file
    type: File
  umap_png:
    outputSource: annotate-concatenate/umap_png
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
      - processed_h5ad_file
      - umap_png

    run: steps/annotate-concatenate.cwl
    label: "Annotates and concatenates h5ad data files in directory"