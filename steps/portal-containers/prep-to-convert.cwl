cwlVersion: v1.0
class: CommandLineTool
label: Renames files for use with portal conversions

hints:
  DockerRequirement:
    dockerPull: hubmap/rna-data-products-python
baseCommand: /opt/prep_to_convert.py

inputs:
  processed_h5ad_file:
    type: File
    doc: The annotated and processed h5ad file
    inputBinding:
      prefix: --processed_h5ad_file

outputs:
  secondary_analysis_dir:
    type: Directory
    outputBinding:
      glob: "processed_dir"