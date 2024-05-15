cwlVersion: v1.0
class: CommandLineTool
label: Annotates each h5ad file with dataset and tissue type, then concatenates

hints:
  DockerRequirement:
    dockerPull: hubmap/rna-data-products
baseCommand: /opt/cell_type_data_products.py

inputs:
  enable_manhole:
    label: "Whether to enable remote debugging via 'manhole'"
    type: boolean?
    inputBinding:
      position: 0

  organ_raw_products:
    type: File[]
    doc: List of paths to raw h5ad files for each organ with cell type annotations
    inputBinding:
      position: 1

outputs:
  raw_h5ad_files:
    type: File[]
    outputBinding:
      glob: "*_raw.h5ad"
    doc: h5ad file containing raw expression data

  processed_h5ad_files:
    type: File[]
    outputBinding:
      glob: "*_processed.h5ad"
    doc: h5ad file containing processed expression data
  
  umap_png:
    type: File[]
    outputBinding:
      glob: "*/*.png"
    doc: leiden cluster umap png