#!/usr/bin/env python3
import json
import os
from argparse import ArgumentParser
from pathlib import Path


def get_uuid(data_product_metadata):
    with open(data_product_metadata, 'r') as json_file:
        metadata = json.load(json_file)
    uuid = metadata["Data Product UUID"]
    return uuid


def upload_to_ec2(umap_png, metadata_json, shiny_cell_dir, uuid, ssh_key):
    os.system(f"scp -i {ssh_key} {umap_png} main_user@ec2-44-213-71-141.compute-1.amazonaws.com:/pipeline_outputs/{uuid}.png")
    os.system(f"scp -i {ssh_key} {metadata_json} main_user@ec2-44-213-71-141.compute-1.amazonaws.com:/pipeline_outputs/{uuid}.json")
    os.system(f"ssh -i {ssh_key} main_user@ec2-44-213-71-141.compute-1.amazonaws.com mkdir {uuid}")
    os.system(f"scp -i {ssh_key} {shiny_cell_dir}/* main_user@ec2-44-213-71-141.compute-1.amazonaws.com:{uuid}")
    os.system(f"ssh -i {ssh_key} main_user@ec2-44-213-71-141.compute-1.amazonaws.com sudo mkdir /srv/shiny-server/{uuid}")
    os.system(f"ssh -i {ssh_key} main_user@ec2-44-213-71-141.compute-1.amazonaws.com sudo mv {uuid}/* /srv/shiny-server/{uuid}")


def main(umap_png, metadata_json, shiny_cell_dir, ssh_key):
    uuid = get_uuid(metadata_json)
    upload_to_ec2(umap_png, metadata_json, shiny_cell_dir, uuid, ssh_key)


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("umap_png", type=Path)
    p.add_argument("data_product_metadata", type=Path)
    p.add_argument("shiny_cell_dir", type=Path)
    p.add_argument("ssh_key", type=Path)
    args = p.parse_args()

    main(args.umap_png, args.data_product_metadata, args.shiny_cell_dir, args.ssh_key)
