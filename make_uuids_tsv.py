#!/usr/bin/env python3
import json
from argparse import ArgumentParser

import pandas as pd
import requests
import yaml

organ_types_yaml_file = "bin/organ_types.yaml"


def extract_donor_metadata(metadata):
    donor_info = {
        "age": None,
        "sex": None,
        "height": None,
        "weight": None,
        "bmi": None,
        "cause_of_death": None,
        "race": None 
    }

    for item in metadata.get('organ_donor_data', []):
        concept = item.get('grouping_concept_preferred_term')
        value = item.get('data_value')
        if concept == "Age":
            donor_info["age"] = value
        elif concept == "Sex":
            donor_info["sex"] = item.get('preferred_term')
        elif concept == "Height":
            donor_info["height"] = value
        elif concept == "Weight":
            donor_info["weight"] = value
        elif concept == "Body mass index":
            donor_info["bmi"] = value
        elif concept == "Cause of death":
            donor_info["cause_of_death"] = item.get('preferred_term')
        elif concept == "Race":
            donor_info["race"] = item.get('preferred_term')

    return donor_info


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
        uuids = []
        hubmap_ids = []
        donor_metadata_list = []

        for hit in hits:
            source = hit["_source"]
            uuids.append(source["uuid"])
            hubmap_ids.append(source["hubmap_id"])
            donor_metadata = source.get("donor", {}).get("metadata", {})
            donor_metadata_list.append(extract_donor_metadata(donor_metadata))

        return uuids, hubmap_ids, donor_metadata_list
    
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
            return uuids, hubmap_ids
    else:
        print("Error:", response.status_code)
        return [], [], []
    

def fill_in_nans(uuids_df: pd.DataFrame) -> pd.DataFrame:
    _donor_cache = {}
    
    for index, row in uuids_df.iterrows():
        if pd.isna(row["age"]) or pd.isna(row["sex"]) or pd.isna(row["race"]):
            uuid = row["uuid"]
            # Check if donor data is already cached
            if uuid not in _donor_cache:
                url = f"https://entity.api.hubmapconsortium.org/datasets/{uuid}/donors"
                response = requests.get(url)
                # Ensure the request is successful
                if response.status_code == 200:
                    donor_data = json.loads(response.content)[0]
                    metadata = donor_data["metadata"]["living_donor_data"]
                    # Populate donor info
                    donor = {}
                    for md in metadata:
                        key = md["grouping_concept_preferred_term"].lower()
                        if key in ["sex", "race", "cause of death"]:
                            donor[key] = md["preferred_term"]
                        elif key == "age":
                            donor[key] = md["data_value"]
                            if md["units"] == "months":
                                donor[key] = str(float(donor[key]) / 12)
                        elif key == "weight":
                            donor[key] = md["data_value"]
                        elif key == "height":
                            donor[key] = md["data_value"]
                        elif key == "body mass index":
                            donor["bmi"] = md["data_value"]
                    _donor_cache[uuid] = donor  # Cache the donor data
                else:
                    print(f"Failed to fetch data for uuid: {uuid}")
                    continue
            else:
                donor = _donor_cache[uuid]
            # Update the DataFrame with donor data
            for key, value in donor.items():
                uuids_df.loc[index, key] = value
    return uuids_df


def main(tissue_type: str):
    organ_dict = yaml.load(open(organ_types_yaml_file), Loader=yaml.BaseLoader)
    for key in organ_dict:
        organ_dict[key] = organ_dict[key]["description"]
    if tissue_type not in organ_dict.values():
        print(f"Tissue {tissue_type} not found ")
        tissue_type = None
    uuids_list, hubmap_ids_list, donor_metadata = get_uuids(organ_dict, tissue_type)
    uuids_df = pd.DataFrame()
    uuids_df["uuid"] = pd.Series(uuids_list, dtype=str)
    uuids_df["hubmap_id"] = pd.Series(hubmap_ids_list, dtype=str)
    donor_metadata_df = pd.DataFrame(donor_metadata)
    result_df = pd.concat([uuids_df, donor_metadata_df], axis=1)
    key_for_tissue = [key for key, value in organ_dict.items() if value == tissue_type]
    if key_for_tissue:
            output_file_name = f"{key_for_tissue[0].lower()}.tsv"
    else:
        output_file_name = "rna.tsv"
    result_df = fill_in_nans(result_df)
    print(result_df)
    result_df.to_csv(output_file_name, sep="\t")


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("tissue_type", type=str, nargs="?", help="Type of tissue (optional)")

    args = p.parse_args()

    main(args.tissue_type)
