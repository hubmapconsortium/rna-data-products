#!/bin/bash
set_access_keys() {
    local access_key_id=$1
    local secret_access_key_id=$2

    aws configure set aws_access_key_id "$access_key_id"
    aws configure set aws_secret_access_key "$secret_access_key"
}
upload_file_to_s3() {
    local local_file=$1
    local remote_path=$2
    file_size=$(stat -c %s "$local_file")
    aws s3 cp "$local_file" "$remote_path" --expected-size "$file_size" --acl public-read
}
upload_files_to_s3(){
    local raw_h5ad=$1
    local processed_h5ad=$2
    local umap_png=$3
    local tissue=$4

    remote_base_path="s3://hubmap-data-products/${tissue:-rna}/"

    upload_file_to_s3 "$raw_h5ad" "$remote_base_path$(basename $raw_h5ad)"
    upload_file_to_s3 "$processed_h5ad" "$remote_base_path$(basename $processed_h5ad)"
    upload_file_to_s3 "$umap_png" "$remote_base_path$(basename $umap_png)"
}

main(){
    local raw_h5ad=$1
    local processed_h5ad=$2
    local umap_png=$3
    local tissue=$4
    local access_key_id=$5
    local secret_access_key=$6

    set_access_keys "$access_key_id" "$secret_access_key"
    upload_files_to_s3 "$raw_h5ad" "$processed_h5ad" "$umap_png" "$tissue"
    echo "Finished" > finished.txt
}

main "$@"
