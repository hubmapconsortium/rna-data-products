#!/usr/bin/env python3
from argparse import ArgumentParser
from pathlib import Path

import anndata
import os
#import s3fs
import zarr
import json
#import awswrangler
import awscli

def set_access_keys(access_key_id: str, secret_access_key: str):
    os.system(f"aws configure set aws_access_key_id \"{access_key_id}\"")
    os.system(f"aws configure set aws_secret_access_key \"{secret_access_key}\"")


def upload_file_to_s3(local_file: Path, uuid: str):
    bucket_path = f"s3://hubmap-data-products/{uuid}/"
    file_size = os.path.getsize(local_file)
    file_name = os.path.basename(local_file)
    os.system(f"aws s3 cp \"{local_file}\" \"{bucket_path}{file_name}\" --expected-size \"{file_size}\"")


def main(zarr_dir: Path, metadata_file: Path, access_key_id: str, secret_access_key: str):
    #stripped_filename = h5ad.removesuffix(".h5ad")
    #adata = anndata.read_h5ad(h5ad)
    #new_adata = anndata.AnnData(obs=adata.obs, obsm=adata.obsm, obsp=adata.obsp)
    #new_adata.write_zarr(f"{stripped_filename}.zarr")
    with open(metadata_file, 'r') as file:
        data = json.load(file)
    uuid = data['Data Product UUID']
    zarr_file = f"{zarr_dir}/secondary_analysis.zarr"
    zarr.consolidate_metadata(zarr_file)
    set_access_keys(access_key_id, secret_access_key)
    upload_file_to_s3(zarr_file, uuid)
    f = open("finished_zarr.txt", "w")
    f.write("zarr uploaded")
    f.close()

if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("zarr_dir", type=Path)
    p.add_argument("metadata_file", type=Path)
    p.add_argument("access_key_id", type=str)
    p.add_argument("secret_access_key", type=str)
    args = p.parse_args()

    main(args.zarr_dir, args.metadata_file, args.access_key_id, args.secret_access_key)
