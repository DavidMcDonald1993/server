import sys
import os.path
sys.path.insert(1, 
    os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))

from pypdb.clients.search.search_client import perform_search
from pypdb.clients.search.search_client import SearchService, ReturnType
from pypdb.clients.search.operators import text_operators

from json.decoder import JSONDecodeError

def get_human_targets():
    print ("getting all targets for Homo sapiens")
    search_service = SearchService.TEXT
    search_operator = text_operators.ExactMatchOperator(value="Homo sapiens",
        attribute="rcsb_entity_source_organism.taxonomy_lineage.name")
    return_type = ReturnType.ENTRY

    human_targets = perform_search(search_service, search_operator, return_type)

    human_targets = {human_target[:4] 
        for human_target in human_targets}
    print ("found", len(human_targets), "human_targets")
    return human_targets

def get_pdb_ids_from_gene_symbol(gene_symbol):
    print ("searching for PDB IDs for gene symbol", gene_symbol)
    search_service = SearchService.TEXT
    search_operator = text_operators.ContainsPhraseOperator(
        value=gene_symbol,
        attribute="struct.title"
        )
    return_type = ReturnType.ENTRY

    try:
        results = perform_search(search_service, search_operator, return_type)
        print ("found", len(results), "results")
        return set(results)

    except JSONDecodeError as e:
        print ("EXCEPTION")
        print (e)
        return set()

if __name__ == "__main__":
    gene_symbol = "CYBA"
    pdb_ids = get_pdb_ids_from_gene_symbol(gene_symbol)
    print (pdb_ids)
