### rna-data-products
A Python and CWL pipeline for concatenating HuBMAP RNA-seq [Salmon] data into data products per organ and one large RNA-seq [Salmon] data product.
## Pipeline steps
* Create a UUIDs TSV file with all UUIDs and HuBMAP IDs of public processed data wanted for the run.
* With the UUIDs TSV, create a data directory of all H5ADs needed for the run.
* Annotate and concatenate a raw data product and a processed data product.
## Requirements
Check the list of python packages in `docker/environment.txt`
## How to run
# Step 1
`python3 make_uuids_tsv.py [tissue_type]`
# Step 2
`python3 make_directory.py /hive/hubmap/data/ [uuids_file]`
# Step 3 
`cwltool pipeline.cwl --[data_directory] --[uuids_file] --[tissue_type]`
