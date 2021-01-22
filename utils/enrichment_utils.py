import os

import sys
import os.path
sys.path.insert(1, 
    os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))

import io

import pandas as pd

from reactomepy.code.rest.service_analysis import AnalysisService
from reactomepy.code.rest.token_mgr import TokenManager

def perform_enrichment_analysis(
    uniprot_id_filename,
    output_csv_filename,
    found_filename,
    not_found_filename,
    pdf_filename=None,
    token=None,
    to_hsa=True,
    resource="TOTAL",
    num_attempts=10,):
    assert isinstance(uniprot_id_filename, str)
    assert os.path.exists(uniprot_id_filename), uniprot_id_filename
   
    print ("performing pathway enrichment")
    ana = AnalysisService()
    print ("generating token using file", uniprot_id_filename)

    for attempt_no in range(num_attempts):

        try:
            token = ana.get_token(uniprot_id_filename,
                token=token,
                to_hsa=to_hsa,)
            break
        except Exception as e:
            print (attempt_no, "REACTOME POST EXCEPTION")
            print (e)
            pass

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

def perform_enrichment_on_uniprot_accs(
    uniprot_confidences, 
    output_dir,
    threshold=500,
    ):

    output_dir = os.path.join(output_dir, "enrichment")
    os.makedirs(output_dir, exist_ok=True)

    print ("perfoming enrichment analysis on uniprot confidences file",
        "to directory", output_dir)

    above_threshold = uniprot_confidences.loc[
        uniprot_confidences["max_confidence"] > threshold]

    unique_uniprots = set(above_threshold["uniprot_ACC"])
    unique_uniprots_filename = os.path.join(output_dir, 
        "unique_uniprot_ACCs.txt")
    print ("writing unique uniprots to", unique_uniprots_filename)
    with open(unique_uniprots_filename, "w") as f:
        f.write("\n".join(unique_uniprots))

    if len(unique_uniprots) > 0:

        # filenames to output enrichment
        output_csv_filename = os.path.join(output_dir, 
            "enrichment.csv")
        found_filename = os.path.join(output_dir,
            "found.txt")
        not_found_filename = os.path.join(output_dir,
            "not_found.txt")
        pdf_filename = os.path.join(output_dir,
            "enrichment_summary.pdf")

        perform_enrichment_analysis(
            unique_uniprots_filename,
            output_csv_filename,
            found_filename,
            not_found_filename,
            pdf_filename)

    return 0

if __name__ == "__main__":
    
    uniprot_id_filename = "user_files/user_id=1/activity_prediction/targets=PARP1_expression_enhancer-thresholds=950-hits-2021-01-22-16:22:15.247338/enrichment/unique_uniprot_ACCs.txt"
    # uniprot_id_filename = "uniprots.txt"
    output_csv_filename = "trash/enrichment.csv"
    found_filename = "trash/found.txt"
    not_found_filename = "trash/not_found.txt"
    pdf_filename = "trash/summary.pdf"

    enrichment, found, not_found =\
        perform_enrichment_analysis(uniprot_id_filename,
        output_csv_filename,
        found_filename, not_found_filename,
        pdf_filename=pdf_filename)

    

