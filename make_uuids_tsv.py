#!/usr/bin/env python3
import json
from argparse import ArgumentParser

import pandas as pd
import requests
import yaml

organ_types_yaml_file = "bin/organ_types.yaml"


def get_uuids(organ_name_mapping: dict, organ: str = None):
    organ_code_mapping = {organ_name_mapping[key]: key for key in organ_name_mapping}
    if organ == None:
        query_payload = {
            "from": 0,
            "size": 10000,
            "query": {
                "bool": {
                    "must": [
                        {"match": {"dataset_type": "RNAseq [Salmon]"}},
                        {"match": {"data_access_level": "public"}},
                    ]
                }
            },
        }
    else:
        query_payload = {
            "from": 0,
            "size": 10000,
            "query": {
                "bool": {
                    "must": [
                        {"match": {"dataset_type": "RNAseq [Salmon]"}},
                        {"match": {"data_access_level": "public"}},
                        {"match": {"origin_samples.organ": organ_code_mapping[organ]}},
                    ]
                }
            },
        }
    query_payload_str = json.dumps(query_payload)
    print(query_payload_str)

    url = "https://search.api.hubmapconsortium.org/v3/search"
    response = requests.post(url, json=query_payload)
    print("Response status code: ", response.status_code)
    print("Response url: ", response.url)
    print("Response reason: ", response.reason)
    if response.status_code == 200:
        data = response.json()
        # Accessing nested dictionaries to retrieve uuid and hubmap_id
        hits = data.get("hits", {}).get("hits", [])
        uuids = [hit["_source"]["uuid"] for hit in hits]
        hubmap_ids = [hit["_source"]["hubmap_id"] for hit in hits]
        data_access_levels = [hit["_source"]["data_access_level"] for hit in hits]

        # Print data_access_level
        print("Data Access Levels:", data_access_levels)
        return uuids, hubmap_ids
    elif response.status_code == 303:
        print("Response body: ", response.text)
        redirection_url = response.text
        print("Redirection URL: ", redirection_url)
        redirection_response = requests.get(redirection_url)
        if redirection_response.status_code == 200:
            data = redirection_response.json()
            hits = data.get("hits", {}).get("hits", [])
            uuids = [hit["_source"]["uuid"] for hit in hits]
            hubmap_ids = [hit["_source"]["hubmap_id"] for hit in hits]
            data_access_levels = [hit["_source"]["data_access_level"] for hit in hits]
            return uuids, hubmap_ids
    else:
        print("Error:", response.status_code)
        return [], []


def main(tissue_type: str):
    organ_dict = yaml.load(open(organ_types_yaml_file), Loader=yaml.BaseLoader)
    for key in organ_dict:
        organ_dict[key] = organ_dict[key]["description"]
    if tissue_type not in organ_dict.values():
        print(f"Tissue {tissue_type} not found ")
        tissue_type = None
    uuids_list, hubmap_ids_list = get_uuids(organ_dict, tissue_type)
    uuids_df = pd.DataFrame()
    uuids_df["uuid"] = pd.Series(uuids_list, dtype=str)
    uuids_df["hubmap_id"] = pd.Series(hubmap_ids_list, dtype=str)
    output_file_name = f"{tissue_type.lower()}.tsv" if tissue_type else "rna.tsv"
    print(uuids_df)
    uuids_df.to_csv(output_file_name, sep="\t")


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("tissue_type", type=str, nargs="?", help="Type of tissue (optional)")

    args = p.parse_args()

    main(args.tissue_type)
