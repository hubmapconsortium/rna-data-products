from os import fspath, walk
from pathlib import Path


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