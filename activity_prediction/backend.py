
import os
import sys
import os.path
sys.path.insert(1, 
    os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))

import pandas as pd

import shutil

from collections import defaultdict

from utils.pass_utils import determine_confidences
from utils.enrichment_utils import perform_enrichment_on_uniprot_accs
from utils.io import download_from_client, get_compound_names, convert_file, write_json
from utils.queries import (get_uniprots_for_targets, get_protein_gene_from_acc, 
    get_drugs_for_uniprots, get_diseases_for_uniprots)
from utils.users import send_file_to_user, determine_identifier
from utils.ppb2_utils import perform_predicton_with_novel_classifier, rescale_predicted_uniprot_confidences
from utils.uniprot_utils import query_uniprot

# def perform_predicton_with_novel_classifier(
#     smiles,
#     model_filename="models/morg3-xgc.pkl.gz",
#     n_proc=6,
#     ):
#     '''
#     Predict from SMILES using novel classifier
#     '''
#     if isinstance(smiles, str): # read from file
#         assert smiles.endswith(".smi")
    
#         # read (and filter smiles)
#         smiles = read_smiles(
#             smiles,
#             filter_valid=True, 
#             return_series=True)

#     assert isinstance(smiles, pd.Series)

#     model = load_model(model_filename)
#     if hasattr(model, "n_proc"):
#         model.set_n_proc(n_proc)
    
#     # make prediction using pretrained model
#     # return as n_targets x n_compounds
#     predictions = model.predict_proba(smiles).T 

#     # id_to_db_id = load_json("id_to_db_id.json")
#     id_to_target_acc = load_json("models/target_ids.json")
    
#     return pd.DataFrame(predictions, 
#         # index=[id_to_db_id[str(i)] for i in range(predictions.shape[1])],
#         index=[id_to_target_acc[str(i)] 
#             for i in range(predictions.shape[0])],
#         columns=smiles.index, 
#     )

# def rescale_predicted_uniprot_confidences(predictions, max_confidence=1000):
#     assert isinstance(predictions, pd.DataFrame)
#     assert predictions.shape[0] == 1683
#     # rescale by max confidence (per compound -- (over all targets)  -- axis 0)
#     return (predictions.divide(predictions.max(axis=0,)) * max_confidence).astype(int)

# def perform_enrichment_on_uniprot_accs(
#     uniprot_confidences, 
#     output_dir,
#     threshold=500,
#     ):

#     output_dir = os.path.join(output_dir, "enrichment")
#     os.makedirs(output_dir, exist_ok=True)

#     print ("perfoming enrichment analysis on uniprot confidences file",
#         "to directory", output_dir)

#     above_threshold = uniprot_confidences.loc[
#         uniprot_confidences["max_confidence"] > threshold]

#     unique_uniprots = set(above_threshold["uniprot_ACC"])
#     unique_uniprots_filename = os.path.join(output_dir, 
#         "unique_uniprot_ACCs.txt")
#     print ("writing unique uniprots to", unique_uniprots_filename)
#     with open(unique_uniprots_filename, "w") as f:
#         f.write("\n".join(unique_uniprots))

#     if len(unique_uniprots) > 0:

#         # filenames to output enrichment
#         output_csv_filename = os.path.join(output_dir, 
#             "enrichment.csv")
#         found_filename = os.path.join(output_dir,
#             "found.txt")
#         not_found_filename = os.path.join(output_dir,
#             "not_found.txt")
#         pdf_filename = os.path.join(output_dir,
#             "enrichment_summary.pdf")

#         perform_enrichment_analysis(
#             unique_uniprots_filename,
#             output_csv_filename,
#             found_filename,
#             not_found_filename,
#             pdf_filename)

#     return 0

def write_actives(confidence_df, threshold, output_dir):
    actives = {compound:
        sorted([(target, confidence) 
            for target, confidence in confidence_df[compound].items()
            if confidence > threshold], key=lambda x: x[1], reverse=True)
            for compound in confidence_df}
    active_targets_filename = os.path.join(output_dir,
        f"actives_(threshold={threshold}).json")
    write_json(actives, active_targets_filename)

def perform_pass_prediction(
    input_file, 
    output_dir, 
    compound_uniprot_confidences,
    threshold=750,
    ):
    '''
    BEGIN PASS
    '''
    assert isinstance(compound_uniprot_confidences, dict)


    '''
    perform pass prediction to directory that will be deleted
    '''
    pass_output_dir = os.path.join(output_dir, "PASS")
    os.makedirs(pass_output_dir, exist_ok=True)

    input_file = convert_file(input_file, 
        desired_format=".sdf", output_dir=pass_output_dir)

    base_name, extension = os.path.splitext(input_file)
    assert extension == ".sdf"

    pass_out_file = base_name + "-PASS-out.sdf"

    out_stream = os.path.join(pass_output_dir, "pass.out")
    err_stream = os.path.join(pass_output_dir, "pass.err")

    cmd = f'''
    PASS2019toSDF.exe {input_file} {pass_out_file}
    '''
    #     > {out_stream} 2> {err_stream}
    # '''
    print ("executing command:", cmd)

    ret = os.system(cmd)
    assert ret == 0
    assert os.path.exists(pass_out_file)

    # DF n_targets * n_compounds
    confidences = determine_confidences(pass_out_file)
    assert isinstance(confidences, dict)

    # delete large SDF file
    print ("deleting pass output directory:", pass_output_dir)
    # os.remove(pass_out_file)
    shutil.rmtree(pass_output_dir)

    # write confidences to file
    # confidence_filename = os.path.join(pass_output_dir, 
    #     "unthresholded_confidences.csv")
    # print ("writing confidences to", confidence_filename)
    # confidences.to_csv(confidence_filename)

    # write jsons for threshold(s)?
    # for threshold in range(500, 1000, 100):
    #     write_actives(confidences, threshold, output_dir=pass_output_dir)

    for compound in confidences:
        compound_output_dir = os.path.join(output_dir,
            compound)
        compound_output_file = os.path.join(compound_output_dir,
            "all_mechanism_toxicity_confidences.tsv")
        compound_target_confidences = confidences[compound]
        targets_sorted_by_confidence = sorted(
            compound_target_confidences, key=compound_target_confidences.get,
            reverse=True)
        with open(compound_output_file, "w") as f:
            f.write("Target\tConfidence\n")
            for target in targets_sorted_by_confidence:
                f.write(f"{target}\t{compound_target_confidences[target]}\n")

    # predicted_uniprot_dir = os.path.join(pass_output_dir, 
    #     "predicted_uniprot_ACCs")
    # os.makedirs(predicted_uniprot_dir, exist_ok=True)

    # compounds_with_no_pass_targets = []

    '''
    associate targets to uniprot ACCs
    '''

    for compound in confidences:
        print ("determining predicted targets for compound", compound,
            "using threshold", threshold)
        compound_confidences = confidences[compound]
        targets = [k for k, v in compound_confidences.items() 
            if v > threshold]

        if len(targets) == 0:
            print ("NO TARGETS PREDICTED FOR COMPOUND", compound, 
                "ABOVE THRESHOLD", threshold)
            # compounds_with_no_pass_targets.append(compound)
            continue

        # get uniprots for targets
        targets_to_uniprot = get_uniprots_for_targets(targets)

        # get proteins / genes from uniprots
        uniprot_confidences = []
        for target, uniprot, protein, gene, association_score in targets_to_uniprot:
            uniprot_confidences.append(
                (target, compound_confidences[target],
                    uniprot, protein, gene, association_score))

        uniprot_confidences = pd.DataFrame(
            sorted(uniprot_confidences, key=lambda x: x[1], reverse=True), # sort by confidence (DESC)
            columns=["target", "target_confidence", "uniprot_ACC", "protein", "gene", "target-acc association_score"])
        uniprot_confidences_filename = os.path.join(
            output_dir, compound,
            f"inferred_uniprot_confidences_thresholded_at_{threshold}.tsv")
        uniprot_confidences.to_csv(uniprot_confidences_filename, sep="\t")

        for _, row in uniprot_confidences.drop_duplicates(
            "uniprot_ACC", keep="first").iterrows():
            uniprot_acc = row["uniprot_ACC"]
            target_confidence = row["target_confidence"]
            # get proteins / genes from uniprot ACC 
            compound_uniprot_confidences[compound].update(
                {uniprot_acc: target_confidence})

    # compounds_with_no_targets_filename = os.path.join(pass_output_dir, 
    #     f"compounds_with_no_predicted_targets_above_{enrichment_threshold}.txt")
    # print ("writing list of missing compounds to", compounds_with_no_targets_filename)
    # with open(compounds_with_no_targets_filename, "w") as f:
    #     f.write("\n".join(compounds_with_no_pass_targets))

    return compound_uniprot_confidences

def perform_ppb2_predict(
    input_file, 
    output_dir, 
    compound_uniprot_confidences,
    model="morg2-nn+nb",
    threshold=750,
    max_confidence=1000):

    '''
    BEGIN PPB2
    '''

    # ppb2_output_dir = os.path.join(output_dir, "PPB2")
    # os.makedirs(ppb2_output_dir, exist_ok=True)

    input_file = convert_file(input_file, 
        desired_format=".smi", output_dir=output_dir)

    # use classifier to generate uniprot confidences
    pred, uniprot_confidences = perform_predicton_with_novel_classifier(
        input_file,
        model=model)

    # pred_filename = os.path.join(
    #     ppb2_output_dir,
    #     "PPB2_uniprot_predictions.csv")
    # print ("writing novel classifier predictions to", 
    #     pred_filename)
    # pred.to_csv(pred_filename)

    # probs_filename = os.path.join(
    #     ppb2_output_dir,
    #     "PPB2_uniprot_probabilities.csv")
    # print ("writing novel classifier probability predictions to", 
    #     probs_filename)
    # uniprot_confidences.to_csv(probs_filename)

    # rescale by compound to range [0, 1000]
    uniprot_confidences = \
        rescale_predicted_uniprot_confidences(
            uniprot_confidences,
        max_confidence=max_confidence)

    # rescaled_probs_filename = os.path.join(ppb2_output_dir,
    #     "PPB2_uniprot_probabilities_rescaled.csv")
    # print ("writing rescaled novel classifier probabilities to", 
    #     rescaled_probs_filename)
    # probs.to_csv(rescaled_probs_filename)

    # write jsons for threshold(s)?
    # for threshold in range(500, 1000, 100):
    #     write_actives(probs, threshold, 
    #         output_dir=ppb2_output_dir)

    accs = list(uniprot_confidences.index)
    acc_to_protein_genes = {
        acc: (protein, genes)
        for acc, protein, genes in get_protein_gene_from_acc(accs)
    }

    for compound in uniprot_confidences: # columns
        print ("determining PPB2 predicted targets for compound", 
            compound,
            "using threshold", threshold)
    
        predicted_compound_uniprot_confidences = uniprot_confidences[compound]

        # compound_confidence_output_dir = os.path.join(ppb2_output_dir, 
        #     "compound_confidences")
        # os.makedirs(compound_confidence_output_dir, exist_ok=True)

        # write acc + protein + gene + confidence to file
        predicted_compound_uniprot_confidences_filename = os.path.join(
            output_dir, compound, 
            f"predicted_uniprot_confidences_thresholded_at_{threshold}.tsv")
        with open(predicted_compound_uniprot_confidences_filename, "w") as f:
            f.write("ACC\tProteinName\tGeneName\tConfidence\n")
            for acc, confidence in predicted_compound_uniprot_confidences.items():
                if confidence < threshold:
                    continue

                # query uniprot to get proteins and genes from ACC
                # if acc not in acc_genes_proteins:
                    # acc_genes_proteins[acc] = query_uniprot(acc=acc)
                protein, genes = acc_to_protein_genes[acc]
                # for protein, gene in protein_genes(acc):
                f.write(f"{acc}\t{protein}\t{genes}\t{confidence}\n")

                # predicted by both models
                if acc in compound_uniprot_confidences[compound]:
                    compound_uniprot_confidences[compound][acc] = max_confidence
                else:
                    compound_uniprot_confidences[compound][acc] = confidence

    return compound_uniprot_confidences

def perform_drug_identificaton(
    compound_uniprot_confidences,
    output_dir):
    '''
    DRUGS
    '''

    # drug_output_dir = os.path.join(output_dir, "drugs")
    # os.makedirs(drug_output_dir, exist_ok=True)

    for compound in compound_uniprot_confidences:
        filename = os.path.join(output_dir, compound,
            "similar_drugs.tsv")

        drugs, cols = get_drugs_for_uniprots(
            list(compound_uniprot_confidences[compound].keys()), as_dict=True)
        
        cols = list(cols)
        idx = cols.index("acc") + 1
        cols.insert(idx, "confidence_score")

        print ("writing drugs for compound", compound, "to file",
            filename)
        with open(filename, "w") as f:
            f.write(f"\t".join(cols) + "\n")
            for drug in drugs:
                acc = drug["acc"]
                confidence_score = compound_uniprot_confidences[compound][acc]
                drug["confidence_score"] = confidence_score
                f.write("\t".join(map(str, (drug[col] for col in cols)))
                    + "\n")

def perform_disease_identification(
    compound_uniprot_confidences,
    output_dir):

    ''' 
    DISEASES
    '''

    # disease_output_dir = os.path.join(output_dir, "diseases")
    # os.makedirs(disease_output_dir, exist_ok=True)

    for compound in compound_uniprot_confidences:
        filename = os.path.join(output_dir, compound,
            "associated_diseases.tsv")

        diseases, cols = get_diseases_for_uniprots(
            list(compound_uniprot_confidences[compound].keys()), 
            as_dict=True)
        print ("writing diseases for compound", compound, "to file",
            filename)

        cols = list(cols)
        idx = cols.index("acc") + 1
        cols.insert(idx, "confidence_score")
        
        with open(filename, "w") as f:
            f.write("\t".join(cols) + "\n")
            for disease in diseases:
                acc = disease["acc"]
                confidence_score = compound_uniprot_confidences[compound][acc]
                disease["confidence_score"] = confidence_score
                f.write("\t".join(map(str, (disease[col] for col in cols)))
                    + "\n")

def activity_predict(
    user,
    input_file, 
    compression="zip",
    root_dir="user_files",
    enrichment_threshold=750,
    max_confidence=1000,
    pass_predict=True,
    ppb2_predict=True,
    model="morg2-nn+nb",
    perform_enrichment=True,
    identify_drugs=True,
    identify_disease=True,
    group_compounds=False,
    ):

    '''
    Perform activity prediction with PASS and PPB2
    '''

    root_dir = os.path.join(root_dir, 
        "user_id={}".format(user.id), "activity_prediction")
    os.makedirs(root_dir, exist_ok=True)

    identifier = determine_identifier(input_file)

    output_dir = os.path.join(root_dir, identifier)
    os.makedirs(output_dir, exist_ok=True)

    # download from client
    input_file = download_from_client(input_file, output_dir)

    # get compound names and write to file
    compound_names = get_compound_names(input_file)
    compound_name_filename = os.path.join(output_dir, "all_compounds.txt")
    print ("writing compound names to", compound_name_filename)
    with open(compound_name_filename, "w") as f:
        f.write("\n".join(compound_names))

    for compound_name in compound_names:
        os.makedirs(os.path.join(output_dir, compound_name),
        exist_ok=True)

    compound_uniprot_confidences = {
        compound_name: dict() for compound_name in compound_names
    }

    if pass_predict:

       compound_uniprot_confidences = perform_pass_prediction(
           input_file, 
           output_dir, 
           compound_uniprot_confidences,
           threshold=enrichment_threshold
       )
    
    if ppb2_predict:
        compound_uniprot_confidences = perform_ppb2_predict(
            input_file,
            output_dir, 
            compound_uniprot_confidences, 
            model=model,
            threshold=enrichment_threshold)
      

    # write joint confidences
    # joint_uniprot_confidences_filename = os.path.join(output_dir, 
        # "all_compound_uniprot_confidences.json")
    # write_json(compound_uniprot_confidences, 
        # joint_uniprot_confidences_filename)

    '''
    write joint predictions to file
    '''

    if group_compounds:
        print ("GROUPING ALL ACCS")
        compound_uniprot_confidences = pd.DataFrame(
            compound_uniprot_confidences)
        max_confidence = compound_uniprot_confidences.max(axis=1)
        compound_uniprot_confidences = {
            "compound_group": max_confidence.to_dict()
        }
        os.makedirs(os.path.join(output_dir, "compound_group"),
            exist_ok=True)
        # joint_uniprot_confidences_filename = os.path.join(output_dir, 
        # "compound_group_uniprot_confidences.json")
        # write_json(compound_uniprot_confidences, 
        #     joint_uniprot_confidences_filename)

    for compound in compound_uniprot_confidences:
        filename = os.path.join(output_dir, compound, 
            "combined_uniprot_confidences.tsv")
        compound_target_confidences = compound_uniprot_confidences[compound]
        with open(filename, "w") as f:
            f.write("ACC\tConfidence\n")
            for target in sorted(compound_target_confidences,
                key=compound_target_confidences.get, reverse=True):
                f.write(f"{target}\t{compound_target_confidences[target]}\n")

    if identify_drugs:
        perform_drug_identificaton(
            compound_uniprot_confidences,
            output_dir)

    if identify_disease:
        perform_disease_identification(
            compound_uniprot_confidences,
            output_dir)

    if perform_enrichment:

        '''
        BEGIN ENRICHMENT ANALYSIS 
        '''

        ret = perform_enrichment_on_uniprot_accs(
            compound_uniprot_confidences,
            output_dir=output_dir, 
            )
    
    # build zip file containing all targets / run settings / run output
    archive_filename = os.path.join(root_dir,
        identifier)
    print ("writing archive to", 
        archive_filename + "." + compression)

    shutil.make_archive(archive_filename, 
        compression, output_dir)

    attachment_filename = f"{archive_filename}.{compression}"
    send_file_to_user(user, attachment_filename, 
        subject="NPAIEngine Activity Prediction and Pathway Enrichment Results")

    return 0

class User:

    def __init__(self, id_):
        self.id = id_

if __name__ == "__main__":
    # name = "david"
    # email = "davemcdonald93@gmail.com"
    # input_file = "/home/david/Desktop/test-PASS-out.sdf"
    user = User(1)
    input_file = "/home/david/Desktop/test.txt"

    ret = activity_predict(user,
        input_file,
        pass_predict=False,
        enrichment=False)

    # # determine_targets(input_file)
    # perform_enrichment_on_PASS_file(input_file, output_dir="/home/david/Desktop", threshold=0)

    # smiles_file = "/home/david/Desktop/pdb_ligands.smi"
    # smiles_file = "coconut_smiles.smi"

    # # read (and filter smiles)
    # smiles = read_smiles(
    #     smiles_file,
    #     filter_valid=True, 
    #     return_series=True)

    # n_compounds = smiles.shape[0]

    # chunksize = 40000
    # n_chunks = n_compounds // chunksize + 1
    # for chunk_no in range(n_chunks):
    #     predictions_filename = f"coconut_uniprot_predictions_chunk_{chunk_no}.csv.gz"
    #     if os.path.exists(predictions_filename):
    #         continue
    #     print ("processing chunk", chunk_no+1)

    #     chunk = smiles[chunk_no*chunksize:(chunk_no+1)*chunksize]

    #     predictions = perform_predicton_with_novel_classifier(chunk)

    #     print ("writing predictions to", predictions_filename)
    #     predictions.to_csv(predictions_filename)