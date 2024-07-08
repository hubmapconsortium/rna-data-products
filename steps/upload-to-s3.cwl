cwlVersion: v1.0
class: CommandLineTool
label: Uploads the annotated and concatenated h5ad files and the umap png

hints:
  DockerRequirement:
    dockerPull: hubmap/rna-data-products-python
baseCommand: /opt/upload_files.py

inputs:
    annotated_raw_h5ad_file:
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

    data_product_metadata
        type: File
        doc: data product metadata json
        inputBinding:
            position: 3
    
    tissue:
        type: string?
        doc: tissue type
        inputBinding: 
            position: 4
    
    access_key_id:
        type: string?
        doc: AWS access key id
        inputBinding:
            position: 5
    
    secret_access_key:
        type: string
        doc: AWS secret access key
        inputBinding:
            position: 6

    ssh_key:
        type: File
        doc: SSH key for the main user of the data products VM
        inputBinding:
            position: 7

outputs: 
    finished_text:
        type: File
        outputBinding:
            glob: "*.txt"