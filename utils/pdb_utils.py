import sys
import os.path
sys.path.insert(1, 
    os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))

from pypdb.clients.search.search_client import perform_search, perform_search_with_graph
from pypdb.clients.search.search_client import SearchService, ReturnType
from pypdb.clients.search.operators import text_operators
from pypdb.clients.search.search_client import QueryGroup, LogicalOperator

from json.decoder import JSONDecodeError

def get_human_targets():
    print ("getting all targets for Homo sapiens")
    search_service = SearchService.TEXT
    search_operator = text_operators.ExactMatchOperator(
        value="Homo sapiens",
        attribute="rcsb_entity_source_organism.taxonomy_lineage.name"
    )
    return_type = ReturnType.ENTRY

    human_targets = perform_search(search_service, search_operator, return_type)

    print ("found", len(human_targets), "human_targets")
    return set(human_targets)

def get_pdb_ids_from_gene_symbol(gene_symbol,
    filter_human=True): # CHANGE FOR UNIPROT DATABASE?
    print ("searching for PDB IDs for gene symbol", gene_symbol)
    search_service = SearchService.TEXT
    # is_human_operator = text_operators.ExactMatchOperator(
        # value="Homo sapiens",
        # attribute="rcsb_entity_source_organism.taxonomy_lineage.name")
    gene_search_operator = text_operators.DefaultOperator(
        value=gene_symbol,
        # attribute="struct.title"
        )
    # QueryGroup associated with being ((Human OR Mus) AND (Under 4 Angstroms))
    # is_human_and_gene_search = QueryGroup(
    #     queries = [is_human_operator, gene_search_operator],
    #     logical_operator = LogicalOperator.AND,
    # )
    return_type = ReturnType.ENTRY

    try:
        results = perform_search(search_service, gene_search_operator, return_type)
        # results = perform_search_with_graph(
            # query_object=is_human_and_gene_search, return_type=return_type)
        if filter_human:
            results = get_human_targets().intersection(results)
        print ("found", len(results), "results")
        return sorted(results)

    except JSONDecodeError as e:
        print ("EXCEPTION")
        print (e)
        return set()

if __name__ == "__main__":
    gene_symbol = "CYBA"
    # gene_symbol = "PARP1"
    pdb_ids = get_pdb_ids_from_gene_symbol(gene_symbol)
    print (pdb_ids)

    # print (get_human_targets())
    # print (human_targets)