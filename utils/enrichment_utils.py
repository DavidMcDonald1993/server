import os

import sys
import os.path
sys.path.insert(1, 
    os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))

import io
from pathlib import Path

import pandas as pd

from reactomepy.code.rest.service_analysis import AnalysisService
from reactomepy.code.rest.token_mgr import TokenManager

def perform_enrichment_analysis(
    uniprot_id_filename,
    output_dir,
    token=None,
    to_hsa=True,
    resource="TOTAL",
    ):
    assert isinstance(uniprot_id_filename, str)
    assert os.path.exists(uniprot_id_filename), uniprot_id_filename


    # filenames to output enrichment
    output_csv_filename = os.path.join(output_dir, 
        "enrichment.csv")
    found_filename = os.path.join(output_dir,
        "found.txt")
    not_found_filename = os.path.join(output_dir,
        "not_found.txt")
    pdf_filename = os.path.join(output_dir,
        "enrichment_summary.pdf")

    print ("outputting csv file to", output_csv_filename)    
    print ("outputting found file to", found_filename)    
    print ("outputting not found file to", not_found_filename)    
    print ("outputting pdf to", pdf_filename)    

    print ("performing pathway enrichment")
    ana = AnalysisService()
    print ("generating token using file", uniprot_id_filename)

    try:
        token = ana.get_token(
            uniprot_id_filename,
            token=token,
            to_hsa=to_hsa,)
    except Exception as e:
        
        print ("Enrichment analysis POST fail", e)
        Path(os.path.join(output_dir, "fail")).touch()
        return

    tm = TokenManager(token)

    if pdf_filename is not None:
        tm.pdf_report(pdf_filename,
            species="Homo sapiens",
            resource=resource)
    
    print ("writing enriched pathways to", output_csv_filename)
    enrichment = tm.csv_pathways(output_csv_filename, 
        resource=resource)
    print ("converting enrichment into dataframe")
    enrichment = pd.read_csv(io.StringIO(enrichment), 
        index_col=0)
    print ("number of enriched pathways:", enrichment.shape[0])

    print ("writing found IDs to",
        found_filename)
    found = tm.csv_found(
        found_filename,
        resource=resource)
    print ("converting cound IDs to dataframe")
    found = pd.read_csv(io.StringIO(found))
    found.index = [str(idx) for idx in found.index]
    print ("number found:", 
        len(set(found["Submitted identifier"])))

    print ("writing not found IDs to",
        not_found_filename)
    not_found = tm.csv_notfound(
        not_found_filename,)
    not_found = not_found.split("\n")[1:]
    print ("number not found:", len(not_found))
    
    return enrichment, found, not_found

def write_uniprots(uniprots, filename):
    print ("writing uniprots to", filename)
    with open(filename, "w") as f:
        f.write("#UNIPROT\n")
        f.write("\n".join(uniprots))


def perform_enrichment_on_uniprot_accs(
    uniprot_confidences, 
    output_dir,
    threshold=500,
    group_compounds=False,
    ):
    assert isinstance(uniprot_confidences, dict)

    output_dir = os.path.join(output_dir, "enrichment")
    os.makedirs(output_dir, exist_ok=True)

    print ("perfoming enrichment analysis on uniprot confidences file",
        "to directory", output_dir)

    if group_compounds:
        uniprot_confidences = pd.DataFrame(uniprot_confidences) # convert to dataframe for ranking
        rank_df = uniprot_confidences.rank(axis=0, ascending=False, na_option="bottom")
        rank_df_filename = os.path.join(output_dir, "target_ranks.csv")
        print ("writing rank df to", rank_df_filename)
        rank_df.to_csv(rank_df_filename)

        mean_rank_df = rank_df.mean(axis=1).sort_values()
        mean_rank_df_filename = os.path.join(output_dir, "mean_ranks.csv")
        print ("writing mean rank df to", mean_rank_df_filename)
        mean_rank_df.to_csv(mean_rank_df_filename)

        n_targets = mean_rank_df.shape[0]
        # just take top 10%
        unique_uniprots = mean_rank_df[:n_targets//10].keys()
        uniprot_filename = os.path.join(output_dir, "uniprot_ACCs.txt")

        write_uniprots(unique_uniprots, uniprot_filename)

        if len(unique_uniprots) > 0:

            perform_enrichment_analysis(
                uniprot_filename,
                output_dir)

    else:
        for compound in uniprot_confidences:
            print ("performing enrichment analysis for compound", compound)

            compound_output_dir = os.path.join(output_dir, 
                str(compound))
            os.makedirs(compound_output_dir, exist_ok=True)

            unique_uniprots = uniprot_confidences[compound].keys()
            uniprot_filename = os.path.join(compound_output_dir, 
                "uniprot_ACCs.txt")
            write_uniprots(unique_uniprots, uniprot_filename)

            if len(unique_uniprots) > 0:

                perform_enrichment_analysis(
                    uniprot_filename,
                    compound_output_dir)
            print ()

    return 0

if __name__ == "__main__":
    
    # uniprot_id_filename = "/home/david/Desktop/enrichment/10-hydroxy aconitine | aconifine /unique_uniprot_ACCs.txt"
    uniprot_id_filename = "unique_uniprot_ACCs.txt"
    # uniprot_id_filename = "uniprots.txt"
    # output_csv_filename = "trash/enrichment.csv"
    # found_filename = "trash/found.txt"
    # not_found_filename = "trash/not_found.txt"
    # pdf_filename = "trash/summary.pdf"

    enrichment, found, not_found =\
        perform_enrichment_analysis(uniprot_id_filename,
        output_dir="trash"
        )
        # output_csv_filename,
        # found_filename, not_found_filename,
        # pdf_filename=pdf_filename)

    

