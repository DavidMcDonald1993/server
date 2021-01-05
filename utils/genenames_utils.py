import sys
import os.path
sys.path.insert(1, 
    os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))

from collections import defaultdict

import httplib2 as http
import json

import re

try:
    from urlparse import urlparse
except ImportError:
    from urllib.parse import urlparse


def prepare_search_terms(search_terms, 
    remove=set(["protein", "agonist", "inhibitor", "antibiotic", "inducer",
        "antagonist", "expression", "enhancer", "treatment", "substrate" ])):
    print ("removing terms", remove, "from search terms")
    search_terms = map(lambda s: 
        re.sub(r"[^a-zA-Z0-9 ]+", " ", s), search_terms) # remove non alpha-numeric characters
    pattern = "({}) *".format("|".join(remove))
    return map(lambda s: re.sub(r"{}".format(pattern), "", s),
        search_terms) # remove words in list

def search_for_hgnc_id(hgnc_id, key="symbol"):
    print ("searching for hgnc_id", hgnc_id)

    headers = {
        'Accept': 'application/json',
    }

    uri = 'http://rest.genenames.org'
    path = "/fetch/hgnc_id/"

    h = http.Http()

    print ("searching url", uri + path + hgnc_id)

    target = urlparse(uri + path + hgnc_id)
    method = 'GET'
    body = ''

    response, content = h.request(
        target.geturl(),
        method,
        body,
        headers)

    if response['status'] == '200':
        # assume that content is a json reply
        # parse content with the json module 
        data = json.loads(content)
        results = data['response']['docs']
        assert len(results) == 1
        result = results[0]
        if key in result:
            return result[key]
        else:
            print ("Key:", key, "not in result")
            return None
    else:
        print( 'Error detected: ' + response['status'])
        return None

def search_for_targets(search_terms, key="symbol", max_hits=100):
    '''
    Use GeneNames REST service 
    returns name, symbol, uniprot
    '''

    headers = {
        'Accept': 'application/json',
    }

    uri = 'http://rest.genenames.org'
    path = "/search/"

    n_search_terms = len(search_terms)
    print ("searching for UNIPROT terms for", 
        n_search_terms, "search term(s):", search_terms)

    search_terms_prepared = prepare_search_terms(search_terms)

    targets = defaultdict(set)

    h = http.Http()

    for term_no, (search_term, search_term_prepared)\
        in enumerate(zip(search_terms, search_terms_prepared)):
        print ("processing term", search_term)
        print ("searching url", uri + path + search_term_prepared)

        target = urlparse(uri + path + search_term_prepared) # handles converting to HTML
        method = 'GET'
        body = ''

        response, content = h.request(
            target.geturl(),
            method,
            body,
            headers)

        if response['status'] == '200':
            # assume that content is a json reply
            # parse content with the json module 
            data = json.loads(content)
            max_score = data["response"]["maxScore"]
            results = data['response']['docs']
            results = list( filter(
                lambda r: r["score"]==max_score, results))
            num_hits = len(results)
            if num_hits > max_hits:
                print ("Too many hits for target", search_term)
                continue # too general
            for i, result in enumerate(results):
                hgnc_id = result["hgnc_id"]
                score = result["score"]

                # fetch info for hgnc id
                data = search_for_hgnc_id(hgnc_id, key=key) # much more data is available
                if data is not None:
                    if not isinstance(data, list):
                        data = [data]
                    for d in data:
                        targets[search_term].add((d, score))
                print ("completed result", i+1, "/",
                    num_hits, "for term no", term_no+1,
                    "/", n_search_terms)
        else:
            print( 'Error detected: ' + response['status'])
        print ("completed term no", term_no+1, 
            "/", n_search_terms)
        print ()

    return targets

def targets_to_gene_name(targets, ):
    if isinstance(targets, list):
        targets = search_for_targets(targets, key="name")
    assert isinstance(targets, dict)
    return {name
        for target in targets
        for name, score in targets[target]
    }

def targets_to_gene_symbol(targets, ):
    if isinstance(targets, list):
        targets = search_for_targets(targets, key="symbol")
    assert isinstance(targets, dict)
    return {symbol
        for target in targets
        for symbol, score in targets[target]
    }

def targets_to_uniprot_ids(targets, ): # loses target->uniprot mapping
    if isinstance(targets, list) or isinstance(targets, set):
        targets = search_for_targets(targets, key="uniprot_ids")
    assert isinstance(targets, dict)
    return {uniprot_id
        for target in targets
        for uniprot_id, score in targets[target]
    }

if __name__ == "__main__":

    # print (targets_to_gene_symbol(
    # ["Alpha 1a adrenoreceptor antagonist"], prepare=True ))
    # print (search_for_hgnc_id("HGNC:3018"))

    targets = ["1,2-alpha-L-fucosidase inhibitor"]

    search_terms = prepare_search_terms(targets)
    print (list(search_terms))

    print (targets_to_gene_name(targets))