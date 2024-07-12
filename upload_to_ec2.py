#!/usr/bin/env python3
import os
from argparse import ArgumentParser
from pathlib import Path


def upload_to_ec2(umap_png, metadata_json, uuid, ssh_key):
    os.system(f"scp -i {ssh_key} {umap_png} main_user@ec2-44-213-71-141.compute-1.amazonaws.com:/pipeline_outputs/{uuid}.png")
    os.system(f"scp -i {ssh_key} {metadata_json} main_user@ec2-44-213-71-141.compute-1.amazonaws.com:/pipeline_outputs/{uuid}.json")


def main(umap_png, metadata_json, uuid, ssh_key):
    upload_to_ec2(umap_png, metadata_json, uuid, ssh_key)


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("umap_png", type=Path)
    p.add_argument("data_product_metadata", type=Path)
    p.add_argument("ssh_key", type=Path)