cwlVersion: v1.0
class: CommandLineTool
label: Adds azimuth annotations to raw data product

requirements:
  DockerRequirement:
    dockerPull: hubmap/rna-data-products-python
baseCommand: /opt/add_azimuth_annotations.py

inputs:
    raw_h5ad_file:
        type: File
        doc: raw h5ad file
        inputBinding:
            prefix: --raw_h5ad_file
            position: 0

    tissue:
        type: string?
        doc: two letter code for tissue type
        inputBinding:
            prefix: --tissue
            position: 1

    metadata_json:
        type: File[]
        doc: version metadata
        inputBinding:
            prefix: --metadata_json
            itemSeparator: " "
            position: 2
    
    annotations_csv:
        type: File[]
        doc: azimuth annotations
        inputBinding:
            prefix: --annotations_csv
            itemSeparator: " "
            position: 3

outputs:
    annotated_raw_h5ad_file:
        type: File
        outputBinding:
            glob: "*_raw.h5ad"
