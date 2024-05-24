cwlVersion: v1.0
class: CommandLineTool
label: Adds azimuth annotations to raw data product

requirements:
  DockerRequirement:
    dockerPull: hubmap/rna-data-products-python
baseCommand: /opt/add_azimuth_annotations.py

inputs:
    metadata_json:
        type: File
        doc: version metadata
        inputBinding:
            position: 0
    
    raw_h5ad_file
        type: File
        doc: raw h5ad file
        inputBinding:
            position: 1
    
    annotations_csv
        type: File
        doc: azimuth annotations
        inputBinding:
            position: 2

    tissue
        type: string?
        doc: two letter code for tissue type
        inputBinding:
            position: 3

outputs:
    annotated_raw_h5ad_file
        type: File
        outputBinding:
            glob: "*_raw.h5ad"
