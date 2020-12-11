import sys
import os.path
sys.path.insert(1, 
    os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))

from utils.io import load_json

def get_categories():
    return {
        "EFFECTS",
        "MECHANISMS",
        'TOXICITY', 
        'ANTITARGETS',
        'METABOLISM', 
        'GENE_EXPRESSION', 
        'TRANSPORTERS'
    }

def get_all_compounds():
    return load_json("compound_ids.json")

def get_all_targets():
    return load_json("target_ids.json")

def get_targets_for_category(category):
    assert category in get_categories()
    return load_json(f"{category}_targets.json",)

def remove_invalid_characters(s):
    return s.replace(".", "")

def parse_pass_spectra(s, mapping=None):
    split = s.split()
    pa = split[0]
    pi = split[1]
    activity = " ".join(split[2:])
    if mapping is not None:
        assert activity in mapping 
        activity = mapping[activity]
    return activity, {"Pa": pa, "Pi": pi}

if __name__ == "__main__":
    for category in get_categories():
        get_targets_for_category(category)

    get_all_targets()