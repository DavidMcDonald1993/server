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
    resource="TOTAL"):
    assert isinstance(uniprot_id_filename, str)
    assert os.path.exists(uniprot_id_filename), uniprot_id_filename
   
    print ("performing pathway enrichment")
    ana = AnalysisService()
    print ("generating token using file", uniprot_id_filename)
    token = ana.get_token(uniprot_id_filename,
        token=token,
        to_hsa=to_hsa,)

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

if __name__ == "__main__":
    
    uniprot_id_filename = "jupyter_notebooks/1q21o3.txt"
    output_csv_filename = "jupyter_notebooks/enrichment.csv"
    found_filename = "jupyter_notebooks/found.txt"
    not_found_filename = "jupyter_notebooks/not_found.txt"
    pdf_filename = "jupyter_notebooks/summary.pdf"

    enrichment, found, not_found =\
        perform_enrichment_analysis(uniprot_id_filename,
        output_csv_filename,
        found_filename, not_found_filename,
        pdf_filename=pdf_filename)

    

