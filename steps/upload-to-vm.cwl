cwlVersion: v1.0
class: CommandLineTool
label: Uploads the annotated and concatenated h5ad files and the umap png

hints:
  DockerRequirement:
    dockerPull: hubmap/rna-data-products-python
baseCommand: /opt/upload_to_vm.py

inputs:
    final_raw_h5mu_file:
        type: File
        doc: The raw h5mu file
        inputBinding:
            position: 0

    processed_h5mu_file:
        type: File
        doc: The processed h5mu file
        inputBinding: 
            position: 1

    umap_png:
        type: File
        doc: PNG of UMAP
        inputBinding: 
            position: 2

    final_data_product_metadata:
        type: File
        doc: data product metadata json
        inputBinding:
            position: 3
    
    zarr_store:
        type: File
        doc: Zipped zarr store of h5ad
        inputBinding:
            position: 4

outputs: 
    finished_text:
        type: File
        outputBinding:
            glob: "*.txt"
