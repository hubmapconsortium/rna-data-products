cwlVersion: v1.0
class: CommandLineTool
label: Generates ShinyCell app

requirements:
  DockerRequirement:
    dockerPull: hubmap/rna-data-products-shinycell

inputs:
  processed_h5ad_file:
    type: File
    doc: The annotated and processed h5ad file
    inputBinding:
      position: 0

  tissue:
    type: string?
    doc: tissue type
    inputBinding:
      position: 1

outputs:
  shinyDir:
    type: Directory
    outputBinding:
      glob: "shinyApps"

baseCommand: [Rscript, /opt/make_shinycell.R]
