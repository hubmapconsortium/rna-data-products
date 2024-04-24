from argparse import ArgumentParser
import awswrangler as wr
import boto3
from pathlib import Path


def upload_file_to_s3(path_to_file, boto_session):
    remote_path = f"s3://hubmap-data-products/{path_to_file.name}"
    with open(path_to_file, "rb") as local_f:
        wr.s3.upload(local_file=local_f, path=remote_path, boto3_session=boto_session)


def upload_files_to_s3(paths_to_files, access_key_id, secret_access_key):
    boto_session = boto3.Session(access_key_id, secret_access_key)
    for path in paths_to_files:
        upload_file_to_s3(path, boto_session)


def main(raw_h5ad, processed_h5ad, umap, access_key_id, secret_access_key):
    files = [raw_h5ad, processed_h5ad, umap]
    upload_files_to_s3(files, access_key_id, secret_access_key)


if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('raw_h5ad_file', type=Path)
    p.add_argument('processed_h5ad_file', type=Path)
    p.add_argument('umap_png', type=Path, nargs='?')
    p.add_argument('access_key_id', type=str)
    p.add_argument('secret_access_key', type=str)


    args = p.parse_args()
    main(args.raw_h5ad_file, args.processed_h5ad_file, args.umap_png, args.access_key_id, args.secret_access_key)