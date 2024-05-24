#!/usr/bin/env python3

"""
Util for copying package & reference metadata from JSON to anndata file in 'uns' slot.
This step could not be bundled with the upstream R script because of an issue with
reticulate. See: https://github.com/rstudio/reticulate/issues/209
"""

from argparse import ArgumentParser
import json
from pathlib import Path
import warnings

import anndata
import pandas as pd

def main(version_metadata: Path, raw_h5ad_file: Path, annotations_csv: Path, tissue_type: str):
    # Load version metadata from JSON
    with open(version_metadata, "rb") as f:
        metadata = json.load(f)

    # Read AnnData file
    ad = anndata.read_h5ad(raw_h5ad_file)

    # Read annotations CSV
    annotations_df = pd.read_csv(annotations_csv)
    annotations_df = annotations_df.set_index('barcodes')

    print(annotations_df.head())

    # Check if annotation was performed
    if metadata["is_annotated"]:
        organ = metadata["azimuth_reference"]["name"]

        if organ in ["lung", "heart", "kidney"]:
            # Load additional metadata
            with open("/data/all_metadata.json", 'r') as j:
                all_data = json.loads(j.read())
            organ_metadata = all_data[organ]

            # Prepare annotation names
            azimuth_annotation_name = "predicted." + organ_metadata["versions"]["azimuth_reference"]["annotation_level"]
            azimuth_label = "azimuth_label"
            azimuth_id = "azimuth_id"
            cl_id = "predicted_CLID"
            standardized_label = "predicted_label"
            match = "cl_match_type"
            score = "prediction_score"
            metadata["annotation_names"] = [azimuth_label, azimuth_id, cl_id, standardized_label, match, score]

            # Check for version mismatch
            if metadata["azimuth_reference"]["version"] != organ_metadata["versions"]["azimuth_reference"]["version"]:
                warnings.warn(
                    f"The Azimuth reference version does not match the "
                    f"Azimuth reference version used to generate the mapping! "
                    f"{metadata['azimuth_reference']['version']} vs "
                    f"{organ_metadata['versions']['azimuth_reference']['version']}"
                )

            # Load mapping data
            mapping_df = pd.read_csv('/data/all_labels.csv')
            organ_annotation = organ + "_" + organ_metadata["versions"]["azimuth_reference"]["annotation_level"]
            mapping_df = mapping_df.loc[mapping_df['Organ_Level'] == organ_annotation]

            # Create mapping dictionary
            keys = mapping_df['A_L'].tolist()
            a_id_map = mapping_df['A_ID'].tolist()
            cl_id_map = mapping_df['CL_ID'].tolist()
            standardized_label_map = mapping_df['Label'].tolist()
            match_map = mapping_df['CL_Match'].tolist()
            mapping_dict = dict(zip(keys, zip(a_id_map, cl_id_map, standardized_label_map, match_map)))

            # Initialize new columns in ad.obs
            ad.obs[azimuth_label] = ""
            ad.obs[azimuth_id] = ""
            ad.obs[cl_id] = ""
            ad.obs[standardized_label] = ""
            ad.obs[match] = ""
            ad.obs[score] = ""

            # Iterate over ad.obs to add annotations
            for idx, row in ad.obs.iterrows():
                if idx in annotations_df.index:
                    annotation = annotations_df.loc[idx, azimuth_annotation_name]
                    ad.obs.at[idx, azimuth_label] = annotation
                    ad.obs.at[idx, score] = annotations_df.at[idx, azimuth_annotation_name + ".score"]

                    # Get the mapping data
                    mapped_values = mapping_dict.get(annotation.strip(), ["other"] * 4)
                    ad.obs.at[idx, azimuth_id] = mapped_values[0]
                    ad.obs.at[idx, cl_id] = mapped_values[1]
                    ad.obs.at[idx, standardized_label] = mapped_values[2]
                    ad.obs.at[idx, match] = mapped_values[3]

            # Add additional metadata
            metadata["CLID"] = {"version": organ_metadata["versions"]["CL_version"]}
            metadata["azimuth_to_CLID_mapping"] = {"version": organ_metadata["versions"]["mapping_version"]}
            for i, key in enumerate(organ_metadata["reviewers"]):
                metadata[f"reviewer{i + 1}"] = key
            metadata["disclaimers"] = {"text": organ_metadata["disclaimer"]}

    # Add metadata to AnnData object and write to file
    ad.uns["annotation_metadata"] = metadata
    ad.write(f"{tissue_type}_raw.h5ad")

if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("version_metadata", type=Path)
    p.add_argument("raw_h5ad_file", type=Path)
    p.add_argument("annotations_csv", type=Path)
    p.add_argument("tissue", type=str)
    args = p.parse_args()

    main(args.version_metadata, args.raw_h5ad_file, args.annotations_csv, args.tissue)
