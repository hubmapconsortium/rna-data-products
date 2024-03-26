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

  uuids_file
    label: "Path to a file containing a list of uuids for the dataset to be indexed"
    type: File

  tissue:
    label: "String description of tissue type"
    type: str?

outputs:

  raw_h5ad_file:
    outputSource: annotate-concatenate/raw_h5ad_file
    type: File
  processed_h5ad_file:
    outputSource: annotate-concatenate/processed_h5ad_file
    type: File


steps:

  - id: annotate-concatenate
    in:
      - id: data_directories
        source: data_directories
      - id: nexus_token
        source: nexus_token
      - id: uuids_file:
        source: uuids_file

    out:
      - raw_h5ad_file
      - processed_h5ad_file

    run: steps/annotate-concatenate.cwl
    label: "Annotates and concatenates h5ad data files in directory"