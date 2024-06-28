#!/usr/bin/env python3

import json
from argparse import ArgumentParser
from collections import defaultdict
from datetime import datetime
from os import fspath, walk, listdir
from pathlib import Path
from typing import Dict, Tuple

import anndata
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
import scipy.sparse
import uuid
import yaml

GENE_MAPPING_DIRECTORIES = [
    Path(__file__).parent.parent / "data",
    Path("/opt/data"),
]


def get_tissue_type(dataset: str) -> str:
    organ_dict = yaml.load(open("/opt/organ_types.yaml"), Loader=yaml.BaseLoader)
    organ_code = requests.get(
        f"https://entity.api.hubmapconsortium.org/dataset/{dataset}/organs/"
    )
    organ_name = organ_dict[organ_code]
    return organ_name.replace(" (Left)", "").replace(" (Right)", "")


def get_inverted_gene_dict():
    inverted_dict = defaultdict(list)
    gene_mapping = read_gene_mapping()
    for ensembl, hugo in gene_mapping.items():
        inverted_dict[hugo].append(ensembl)
    return inverted_dict


def find_files(directory, patterns):
    for dirpath_str, dirnames, filenames in walk(directory):
        dirpath = Path(dirpath_str)
        for filename in filenames:
            filepath = dirpath / filename
            for pattern in patterns:
                if filepath.match(pattern):
                    return filepath


def find_file_pairs(directory):
    filtered_patterns = ["cluster_marker_genes.h5ad", "secondary_analysis.h5ad"]
    unfiltered_patterns = ["out.h5ad", "expr.h5ad"]
    filtered_file = find_files(directory, filtered_patterns)
    unfiltered_file = find_files(directory, unfiltered_patterns)
    return filtered_file, unfiltered_file


def annotate_file(
    filtered_file: Path, unfiltered_file: Path, tissue_type: str, uuids_df: pd.DataFrame
) -> Tuple[anndata.AnnData, anndata.AnnData]:
    # Get the directory
    data_set_dir = fspath(unfiltered_file.parent.stem)
    # And the tissue type
    tissue_type = tissue_type if tissue_type else get_tissue_type(data_set_dir)
    hubmap_id = uuids_df.loc[uuids_df["uuid"] == data_set_dir, "hubmap_id"].values[0]
    age = uuids_df.loc[uuids_df["uuid"] == data_set_dir, "age"].values[0]
    sex = uuids_df.loc[uuids_df["uuid"] == data_set_dir, "sex"].values[0]
    height = uuids_df.loc[uuids_df["uuid"] == data_set_dir, "height"].values[0]
    weight = uuids_df.loc[uuids_df["uuid"] == data_set_dir, "weight"].values[0]
    bmi = uuids_df.loc[uuids_df["uuid"] == data_set_dir, "bmi"].values[0]
    cause_of_death = uuids_df.loc[uuids_df["uuid"] == data_set_dir, "cause_of_death"].values[0]
    race = uuids_df.loc[uuids_df["uuid"] == data_set_dir, "race"].values[0]
    filtered_adata = anndata.read_h5ad(filtered_file)
    unfiltered_adata = anndata.read_h5ad(unfiltered_file)
    unfiltered_copy = unfiltered_adata.copy()
    unfiltered_copy.obs["barcode"] = unfiltered_adata.obs.index
    unfiltered_copy.obs["dataset"] = data_set_dir
    unfiltered_copy.obs["hubmap_id"] = hubmap_id
    unfiltered_copy.obs["organ"] = tissue_type
    unfiltered_copy.obs["modality"] = "rna"
    unfiltered_copy.obs["age"] = age
    unfiltered_copy.obs["sex"] = sex
    unfiltered_copy.obs["height"] = height
    unfiltered_copy.obs["weight"] = weight
    unfiltered_copy.obs["bmi"] = bmi
    unfiltered_copy.obs["cause_of_death"] = cause_of_death
    unfiltered_copy.obs["race"] = race
    unfiltered_copy.uns["annotation_metadata"] = (
        filtered_adata.uns["annotation_metadata"]
        if "annotation_metadata" in filtered_adata.uns.keys()
        else {"is_annotated": False}
    )
    
    cell_ids_list = [
        "-".join([data_set_dir, barcode]) for barcode in unfiltered_copy.obs["barcode"]
    ]
    unfiltered_copy.obs["cell_id"] = pd.Series(
        cell_ids_list, index=unfiltered_copy.obs.index, dtype=str
    )
    unfiltered_copy.obs.set_index("cell_id", drop=True, inplace=True)
    unfiltered_copy = map_gene_ids(unfiltered_copy)
    return unfiltered_copy


def read_gene_mapping() -> Dict[str, str]:
    """
    Try to find the Ensembl to HUGO symbol mapping, with paths suitable
    for running this script inside and outside a Docker container.
    :return:
    """
    for directory in GENE_MAPPING_DIRECTORIES:
        mapping_file = directory / "ensembl_to_symbol.json"
        if mapping_file.is_file():
            with open(mapping_file) as f:
                return json.load(f)
    message_pieces = ["Couldn't find Ensembl → HUGO mapping file. Tried:"]
    message_pieces.extend(f"\t{path}" for path in GENE_MAPPING_DIRECTORIES)
    raise ValueError("\n".join(message_pieces))


def map_gene_ids(adata):
    obsm = adata.obsm
    uns = adata.uns
    gene_mapping = read_gene_mapping()
    has_hugo_symbol = [gene in gene_mapping for gene in adata.var.index]
    # adata = adata[:, has_hugo_symbol]
    temp_df = pd.DataFrame(
        adata.X.todense(), index=adata.obs.index, columns=adata.var.index
    )
    aggregated = temp_df.groupby(level=0, axis=1).sum()
    adata = anndata.AnnData(aggregated, obs=adata.obs)
    adata.var["hugo_symbol"] = [
        gene_mapping.get(var, np.nan) for var in adata.var.index
    ]
    adata.obsm = obsm
    adata.uns = uns
    # This introduces duplicate gene names, use Pandas for aggregation
    # since anndata doesn't have that functionality
    adata.X = scipy.sparse.csr_matrix(adata.X)
    adata.var_names_make_unique()
    return adata


def split_and_save(adata, base_file_name, max_cells: int = 30000):
    num_cells = adata.n_obs
    num_splits = (num_cells // max_cells) + 1
    for i in range(num_splits):
        start_idx = i * max_cells
        end_idx = min((i + 1) * max_cells, num_cells)
        adata_sub = adata[start_idx:end_idx].copy()
        part_file_name = f"{base_file_name}_part_{i+1}.h5ad"
        adata_sub.write(part_file_name)


def create_json(tissue, data_product_uuid, creation_time, uuids, hbmids):
    bucket_url = f"https://hubmap-data-products.s3.amazonaws.com/{data_product_uuid}/"
    metadata = {
        "Data Product UUID": data_product_uuid,
        "Tissue": tissue,
        "Raw URL": bucket_url + f"{tissue}_raw.h5ad",
        "Processed URL": bucket_url + f"{tissue}_processed.h5ad",
        "Creation Time": creation_time,
        "Dataset UUIDs": uuids,
        "Dataset HBMIDs": hbmids
    }
    print("Writing metadata json")
    with open(f"{tissue}_metadata.json", "w") as outfile:
        json.dump(metadata, outfile)


def main(data_directory: Path, uuids_file: Path, tissue: str = None):
    raw_output_file_name = f"{tissue}_raw" if tissue else "rna_raw"
    uuids_df = pd.read_csv(uuids_file, sep="\t", dtype=str)
    uuids_list = uuids_df["uuid"]
    hbmids_list = uuids_df["hubmap_id"]
    directories = [data_directory / Path(uuid) for uuid in uuids_df["uuid"]]
    # Load files
    file_pairs = [find_file_pairs(directory) for directory in directories if len(listdir(directory))>1]
    print("Annotating objects")
    adatas = [
        annotate_file(file_pair[0], file_pair[1], tissue, uuids_df)
        for file_pair in file_pairs
    ]
    annotation_metadata = {
        adata.obs.dataset.iloc[0]: adata.uns["annotation_metadata"] for adata in adatas
    }
    saved_var = adatas[0].var
    print("Concatenating objects")
    adata = anndata.concat(adatas, join="outer")
    adata.uns["annotation_metadata"] = annotation_metadata
    creation_time = str(datetime.now())
    adata.uns["creation_date_time"] = creation_time
    adata.uns["datasets"] = hbmids_list
    data_product_uuid = uuid.uuid4()
    adata.uns["uuid"] = uuid.uuid4()
    adata.var = saved_var
    print(f"Writing {raw_output_file_name}")
    adata.write(f"{raw_output_file_name}.h5ad")
    print(f"Splitting and writing {raw_output_file_name} into multiple files")
    split_and_save(adata, raw_output_file_name)
    create_json(tissue, data_product_uuid, creation_time, uuids_list, hbmids_list)


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("data_directory", type=Path)
    p.add_argument("uuids_file", type=Path)
    p.add_argument("tissue", type=str, nargs="?")
    p.add_argument("--enable_manhole", action="store_true")

    args = p.parse_args()

    if args.enable_manhole:
        import manhole

        manhole.install(activate_on="USR1")

    main(args.data_directory, args.uuids_file, args.tissue)