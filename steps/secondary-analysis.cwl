cwlVersion: v1.0
class: CommandLineTool
label: Perform secondary analysis on raw data product

hints:
  DockerRequirement:
    dockerPull: hubmap/rna-data-products-python
baseCommand: /opt/azimuth.R

inputs: 
    raw_h5ad_file
        type: File
        doc: The raw h5ad file
        inputBinding:
            position: 0
    
    annotated_csv
        type: File
        doc: csv with azimuth annotations
        inputBinding:
            position: 1

    tissue:
        type: str?
        doc: optional tissue type
        inputBinding
            position: 2

outputs:
    annotated_raw_h5ad_file:
        type: File
        outputBinding:
            glob: "*_raw.h5ad"
        doc: raw h5ad file with azimuth annotations
    
    processed_h5ad_file:
        type: File
        outputBinding:
            glob: "*_processed.h5ad"
        doc: h5ad file with secondary analysis processing
    
    umap_png
        type: File
        outputBinding:
            glob: *.png
        doc: umap png