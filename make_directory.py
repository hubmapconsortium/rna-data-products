#!/usr/bin/env python3
import json
from argparse import ArgumentParser
from os import fspath, walk
from pathlib import Path
from subprocess import check_call
import pandas as pd


def find_files(directory, patterns):
    for dirpath_str, dirnames, filenames in walk(directory):
        dirpath = Path(dirpath_str)
        for filename in filenames:
            filepath = dirpath / filename
            for pattern in patterns:
                if filepath.match(pattern):
                    return filepath


def find_file_pairs(directory):
    filtered_patterns = ['cluster_marker_genes.h5ad', 'secondary_analysis.h5ad']
    unfiltered_patterns = ['out.h5ad', 'expr.h5ad']
    filtered_file = find_files(directory, filtered_patterns)
    unfiltered_file = find_files(directory, unfiltered_patterns)
    return filtered_file, unfiltered_file


def get_input_directory(data_directory, uuid):
    public_directory = data_directory / 'public' / uuid
    if public_directory.exists():
        return public_directory
    else:
        consortium_directory = data_directory / 'consortium'
        if consortium_directory.exists():
            for subdir in consortium_directory.iterdir():
                consortium_subdir = subdir / uuid
                if consortium_subdir.exists():
                    return consortium_subdir


def main(data_directory: Path, uuids_file: Path):
    uuids = pd.read_csv(uuids_file, sep='\t')[["uuid"], ["hubmap_id"]]
    uuids = uuids.dropna()
    h5ads_base_directory = Path("h5ads")
    h5ads_base_directory.mkdir(exist_ok=True)  # Create h5ads directory if it doesn't exist
    for uuid in uuids["uuid"]:
        h5ads_directory = h5ads_base_directory / uuid
        h5ads_directory.mkdir(parents=True, exist_ok=True)  # Create UUID-specific directory
        input_directory = get_input_directory(data_directory, uuid)
        input_files = find_file_pairs(input_directory)
        if input_files == (None, None):
            print("No input files in: ", input_directory)
            continue
        print("Input directory:", input_directory)
        print("Input files:",input_files)
        for input_file in input_files:
            check_call(f"cp {fspath(input_file)} {h5ads_directory}/{input_file.name}", shell=True)


if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('data_directory', type=Path)
    p.add_argument('uuids_file', type=Path)

    args = p.parse_args()

    main(args.data_directory, args.uuids_file)
