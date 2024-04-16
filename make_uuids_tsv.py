#!/usr/bin/env python3
import yaml
import requests
from argparse import ArgumentParser
import pandas as pd

organ_types_yaml_file = 'bin/organ_types.yaml'

def get_uuids(organ_name_mapping:dict, organ:str=None):
  organ_code_mapping = {organ_name_mapping[key]:key for key in organ_name_mapping}
  url = f'https://search.api.hubmapconsortium.org/v3/param-search/datasets?dataset_type=RNAseq%20%5BSalmon%5D&origin_samples.organ={organ_code_mapping[organ]}' if organ else \
    'https://search.api.hubmapconsortium.org/v3/param-search/datasets?dataset_type=RNAseq%20%5BSalmon%5D'
  uuids = [item['uuid'] for item in requests.get(url).json()]
  hubmap_ids = [item['hubmap_id'] for item in requests.get(url).json()]
  return uuids, hubmap_ids

def main(tissue_type: str):
    organ_dict = yaml.load(open(organ_types_yaml_file), Loader=yaml.BaseLoader)
    for key in organ_dict:
        organ_dict[key] = organ_dict[key]['description']

    if tissue_type not in organ_dict.values():
        print(f"Tissue {tissue_type} not found ")
        tissue_type = None

    uuids_list, hubmap_ids_list = get_uuids(organ_dict, tissue_type)
    uuids_df = pd.DataFrame()
    uuids_df['uuid'] = pd.Series(uuids_list)
    uuids_df['hubmap_id'] = pd.Series(hubmap_ids_list)
    output_file_name = f"{tissue_type.lower()}.tsv" if tissue_type else "rna.tsv"
    uuids_df.to_csv(output_file_name, sep='\t')

if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('tissue_type', type=str)

    args = p.parse_args()

    main(args.tissue_type)
