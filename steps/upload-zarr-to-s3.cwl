cwlVersion: v1.0
class: CommandLineTool
label: Uploads the zarr files to s3

hints:
  DockerRequirement:
    dockerPull: hubmap/rna-data-products-python
baseCommand: /opt/upload_zarr_to_s3.py

inputs:
    zarr_dir:
        type: Directory
        doc: The zarr files
        inputBinding:
            position: 0

    metadata_file:
        type: File
        doc: data product metadata json
        inputBinding:
            position: 1
    
    access_key_id:
        type: string
        doc: AWS access key id
        inputBinding:
            position: 2
    
    secret_access_key:
        type: string
        doc: AWS secret access key
        inputBinding:
            position: 3

outputs: 
    finished_text:
        type: File
        outputBinding:
            glob: "*.txt"
