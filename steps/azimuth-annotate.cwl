cwlVersion: v1.0
class: CommandLineTool
label: Adds azimuth annotations to raw data product

hints:
  DockerRequirement:
    dockerPull: hubmap/rna-data-products-r 
baseCommand: /opt/azimuth.R

inputs:
    raw_h5ad_file:
        type: File
        doc: The raw h5ad file
        inputBinding:
            position: 0
    tissue:
        type: str?
        doc: tissue type
        inputBinding:
            position: 1

outputs: 
    annotations_csv
        type: File?
        outputBinding:
            glob: "*.csv"