cwlVersion: v1.0
class: CommandLineTool
label: Uploads the annotated and concatenated h5ad files and the umap png

hints:
  DockerRequirement:
    dockerPull: hubmap/rna-data-products
baseCommand: /opt/upload_to_s3.py

inputs:
    raw_h5ad_file:
        type: File
        doc: The raw h5ad file
        inputBinding:
            position: 0

    processed_h5ad_file:
        type: File
        doc: The processed h5ad file
        inputBinding: 
            position: 1

    umap_png:
        type: File
        doc: PNG of UMAP
        inputBinding: 
            position: 2
    
    access_key_id:
        type: string
        doc: AWS access key id
        inputBinding: 
            position: 3
    
    secret_access_key:
        type: string
        doc: AWS secret access key
        inputBinding: 
            position: 4