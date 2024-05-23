cwlVersion: v1.0
class: CommandLineTool
label: Adds azimuth annotations to raw data product

requirements:
  DockerRequirement:
    dockerPull: hubmap/rna-data-products-azimuth

inputs:
    raw_h5ad_file:
        type: File
        doc: The raw h5ad file
        inputBinding:
            position: 0

    tissue:
        type: string?
        doc: tissue type
        inputBinding:
            position: 1

outputs: 
    annotations_csv:
        type: File?
        outputBinding:
            glob: "*.csv"
    metadata_json:
        type: File
        outputBinding:
            glob: "*.json"

baseCommand: [Rscript, /opt/azimuth.R]