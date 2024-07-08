from argparse import ArgumentParser
import json
import os


def set_access_keys(access_key_id, secret_access_key):
    os.system(f"aws configure set aws_access_key_id \"{access_key_id}\"")
    os.system(f"aws configure set aws_secret_access_key \"{secret_access_key}\"")


def upload_file_to_s3(local_file, remote_path):
    file_size = os.path.getsize(local_file)
    os.system(f"aws s3 cp \"{local_file}\" \"{remote_path}\" --expected-size \"{file_size}\"")


def upload_files_to_s3(local_paths, remote_paths):

    for local_file, remote_file in local_paths, remote_paths:
        upload_file_to_s3(local_file, remote_file)


def upload_to_ec2(umap_png, metadata_json, uuid, ssh_key):
    os.system(f"scp -i {ssh_key} {umap_png} main_user@ec2-44-213-71-141.compute-1.amazonaws.com:/pipeline_outputs/{uuid}.png")
    os.system(f"scp -i {ssh_key} {metadata_json} main_user@ec2-44-213-71-141.compute-1.amazonaws.com:/pipeline_outputs/{uuid}.json")


def get_uuid(metadata_json):
    with open(metadata_json) as json_file:
        metadata = json.load(json_file)
    uuid = metadata["Data Product UUID"]
    return uuid


def main(raw_h5ad, processed_h5ad, umap_png, data_product_metadata, tissue, access_key_id, secret_access_key, ssh_key):
    set_access_keys(access_key_id, secret_access_key)
    uuid = get_uuid(data_product_metadata)
    bucket_path = f"s3://hubmap-data-products/{uuid}/"
    raw_h5ad_path = bucket_path+f"{tissue}_raw.h5ad"
    processed_h5ad_path = bucket_path+f"{tissue}_processed.h5ad"
    umap_path = bucket_path+f"{tissue}.png"
    remote_file_paths = [raw_h5ad_path, processed_h5ad_path, umap_path]
    local_file_paths = [raw_h5ad, processed_h5ad, umap_png]
    upload_files_to_s3(local_file_paths, remote_file_paths)
    upload_to_ec2(umap_png, data_product_metadata, uuid, ssh_key)
    f = open("finished.txt", "w")
    f.write("cwl wants an output file for this step")
    f.close()


# def main(umap, metadata, access_key_id, secret_access_key, ssh_key, ec2_user):
#     set_access_keys(access_key_id, secret_access_key)
#     uuid = get_uuid(metadata)
#     upload_to_ec2(umap, metadata, uuid, ssh_key, ec2_user)


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("annotated_raw_h5ad_file", type=str)
    p.add_argument("processed_h5ad_file", type=str)
    p.add_argument("umap_png", type=str)
    p.add_argument("data_product_metadata", type=str)
    p.add_argument("tissue", type=str)
    p.add_argument("access_key_id", type=str)
    p.add_argument("secret_access_key", type=str)
    p.add_argument("ssh_key")
    args= p.parse_args()

    main(args.annotated_raw_h5ad_file, args.processed_h5ad_file, args.umap_png, args.data_product_metadata, args.tissue, args.access_key_id, args.secret_access_key, args.ssh_key)