import sys
import os.path
sys.path.insert(1, 
    os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))

import requests
import json

import urllib.parse
import urllib.request


def query_uniprot(acc=None, gene=None, size=1000):

    if acc.startswith("CHEMBL"):
        print ("invalid ACC", acc) #TODO
        return []

    if acc is not None:
        request_url = f"https://www.ebi.ac.uk/proteins/api/proteins?offset=0&size={size}&accession={acc}"
    else:
        assert gene is not None
        request_url = f"https://www.ebi.ac.uk/proteins/api/proteins?offset=0&size={size}&gene={gene}"

    print ("Querying Uniprot using URL", request_url)

    r = requests.get(request_url, headers={ "Accept" : "application/json"})

    if not r.ok:
        r.raise_for_status()
        sys.exit()

    response_body = r.text
    response_body = json.loads(response_body)
    # if len(response_body) == 1:
        # response_dict = response_body[0]
    hits = []
    for response_dict in response_body:
        assert isinstance(response_dict, dict)
        '''
        dict_keys(['accession', 'id', 'proteinExistence', 'info', 'organism', 
        'secondaryAccession', 'protein', 'gene', 'comments', 'features', 'dbReferences', '
            keywords', 'references', 'sequence'])
        '''
        # print (response_dict["protein"].keys())
        if "protein" in response_dict:
            protein_name = response_dict["protein"]
            if "recommendedName" in protein_name:
                protein_name = protein_name["recommendedName"]
            if "submittedName" in protein_name:
                protein_name = protein_name["submittedName"][0]
            if "fullName" in protein_name:
                protein_name = protein_name["fullName"]["value"]
            else:
                assert False, protein_name
        else:
            protein_name = ""
        if "gene" in response_dict:
            genes = response_dict["gene"]
            assert isinstance(genes, list)
            if "name" in genes[0]:
                genes = [gene["name"]["value"] for gene in genes]
            elif "orfNames" in genes[0]:
                genes = [gene["orfNames"][0]["value"] for gene in genes]
        else:
            genes = [""]
        hits.extend([(protein_name, gene) for gene in genes])
    return hits


def map_to_uniprot(query, map_from="CHEMBL_ID", map_to="ACC", return_format="tab"):
    assert isinstance(query, list)

    print ("mapping query", query, "to UNIPROT")

    params = {
        "from": map_from,
        "to": map_to,
        'format': return_format,
        'query': " ".join(query),
    }

    url = 'https://www.uniprot.org/uploadlists/'

    data = urllib.parse.urlencode(params)
    data = data.encode('utf-8')
    req = urllib.request.Request(url, data)
    with urllib.request.urlopen(req) as f:
        response = f.readlines()
    if return_format == "tab":
        delimiter = "\t"
    else:
        raise NotImplementedError
    response = map(lambda s: s.decode("utf-8").rstrip().split(delimiter), response)
    next(response) # drop FROM TO line
    return {chembl: uniprot 
        for chembl, uniprot in response }

if __name__ == "__main__":

    # uniprots = []
    # with open("unique_uniprot_ACCs.txt", "r") as f:
    #     for line in f.readlines():
    #         uniprots.append(line.rstrip())
    # print (uniprots)
    acc = "C4YTQ8"
    # for acc in uniprots:
    hits = query_uniprot(acc=acc)
    print (acc, hits)
    # print (hits)
    # query_uniprot(size=-1, gene="CHEMBL2169726")#gene="PARP1", size=-1)
    # map_to_uniprot(["Q25856"])

    # from utils.io import load_json

    # chembl_targets = load_json("models/id_to_chembl.json")

    # mapping = map_to_uniprot(list(chembl_targets.values()))

    # id_to_uniprot = {
    #     i: (mapping[chembl] if chembl in mapping else chembl)
    #     for i, chembl in chembl_targets.items()
    # }

    # from utils.io import write_json
    # write_json(id_to_uniprot, "models/id_to_uniprot.json")
    