import json 

def load_json(json_filename):
    with open(json_filename, "r") as f:
        return json.load(f)