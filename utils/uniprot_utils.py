import sys
import os.path
sys.path.insert(1, 
    os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))

import requests
import json

def query_uniprot(acc=None, gene=None, size=1000):

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
    print (len(response_body))
    if len(response_body) == 1:
        response_dict = response_body[0]
    assert isinstance(response_dict, dict)
    '''
    dict_keys(['accession', 'id', 'proteinExistence', 'info', 'organism', 
    'secondaryAccession', 'protein', 'gene', 'comments', 'features', 'dbReferences', '
        keywords', 'references', 'sequence'])
    '''
    print (response_dict.keys())

if __name__ == "__main__":
    acc = "Q9NRK6"
    query_uniprot(gene="PARP1", size=-1)