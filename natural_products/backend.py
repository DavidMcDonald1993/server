import os
import sys
import os.path
sys.path.insert(1, 
    os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))

import re

import numpy as np
import pandas as pd

from datetime import datetime

from rdkit import Chem
from rdkit.Chem.Draw import MolToImage

import urllib.parse as urlparse

from utils.io import load_json, write_json, write_smiles
from utils.mongodb_utils import connect_to_mongodb
from utils.mysql_utils import mysql_query, mysql_create_table, mysql_insert_many, connect_to_mysqldb, sanitise_names
from utils.pass_utils import get_categories, get_all_targets, get_targets_for_category, get_all_compounds
from utils.enrichment_utils import perform_enrichment_analysis

def query_target_hits(
    targets, 
    thresholds,
    filter_pa_pi=True,
    limit=None
    ):
    assert filter_pa_pi
    assert isinstance(targets, list)
    if not isinstance(thresholds, list):
        assert isinstance(thresholds, int), thresholds
        # assert thresholds in (900, )
        thresholds = [thresholds]
    num_targets = len(targets)
    if len(thresholds) < num_targets:
        thresholds = thresholds * num_targets

    print ("querying PASS database for compounds",
        "that hit targets", targets,
         "with Pa greater than or equal to",
        "thresholds", thresholds)

    target_names = sanitise_names(targets)
    columns = ", ".join((
        # `{target}_activity`.Pa AS `{target}-Pa`, `{target}_activity`.Pi AS `{target}-Pi`, 
        f'''
            `{target}_activity`.Pa-`{target}_activity`.Pi AS `{target}-Confidence Score`
        '''
        for target in target_names
    ))
    tables = "\n".join((
            f'''
            INNER JOIN activities AS `{target}_activity` 
                ON (c.compound_id=`{target}_activity`.compound_id)
            INNER JOIN targets AS `{target}_target` 
                ON (`{target}_activity`.target_id=`{target}_target`.target_id)
            '''
        for target in target_names
    ))
    conditions = "\n".join((
        f'''
        AND `{target_name}_activity`.above_{threshold}=(1) 
        AND `{target_name}_target`.target_name="{target}"
        '''
        for target_name, target, threshold in zip(target_names[1:], targets[1:], thresholds[1:])      
    ))
    query = f'''
        SELECT c.compound_id, c.coconut_id AS `ID`, c.image_path AS `Image`, c.name AS `Molecule Name`, 
            c.formula AS `Molecular Formula`, c.clean_smiles AS `SMILES`,
        {columns}
        FROM compounds AS c
        {tables}
        WHERE `{target_names[0]}_activity`.above_{thresholds[0]}=(1)
        AND `{target_names[0]}_target`.target_name="{targets[0]}"
        {conditions}
        {f"LIMIT {limit}" if limit is not None else ""}
    '''

    hits, cols = mysql_query(query, return_cols=True)

    records = [
        record[1:] for record in hits
    ]

    return records, cols[1:]

def query_pathway_hits(
    pathways,
    organisms, 
    threshold=0, 
    filter_pa_pi=True, 
    limit=None,
    existing_conn=None):
    assert filter_pa_pi

    if not isinstance(pathways, list):
        assert isinstance(pathways, str)
        pathways = [pathways]
    if not isinstance(organisms, list):
        assert isinstance(organisms, str)
        organisms = [organisms]

    pathway_names = sanitise_names(
        [   f"{pathway}_{organism}"
            for pathway, organism in zip(pathways, organisms)]
    )

    columns = ", ".join((
        f'''
        GROUP_CONCAT(DISTINCT(`{pathway}_target`.target_name) SEPARATOR '-') AS `{pathway} Target Names`, 
        COUNT(DISTINCT(`{pathway}_target`.target_name)) AS `{pathway} Number Targets`, 
        GROUP_CONCAT(DISTINCT(`{pathway}_uniprot`.acc) SEPARATOR '-') AS `{pathway} Uniprot ACCs`, 
        COUNT(DISTINCT(`{pathway}_uniprot`.acc)) AS `{pathway} Number Uniprot ACCs`, 
        `{pathway}_counts`.counts AS `{pathway} Total Uniprot ACCs`,
        CAST(COUNT(DISTINCT(`{pathway}_uniprot`.acc)) / `{pathway}_counts`.counts AS CHAR)
            AS `{pathway} Uniprot Coverage`,   
        `{pathway}_pathway`.pathway_name AS `{pathway} Name`,
        `{pathway}_pathway`.organism AS `{pathway} Organism`,
        `{pathway}_pathway`.pathway_url AS `{pathway} URL`
        '''
        for pathway in pathway_names
    ))

    tables = "\n".join((
        f'''
            INNER JOIN activities AS `{pathway}_activity` 
                ON (c.compound_id=`{pathway}_activity`.compound_id)
            INNER JOIN targets AS `{pathway}_target` 
                ON (`{pathway}_activity`.target_id=`{pathway}_target`.target_id)
            INNER JOIN targets_to_uniprot AS `{pathway}_targets_to_uniprot` 
                ON (`{pathway}_activity`.target_id=`{pathway}_targets_to_uniprot`.target_id)
            INNER JOIN uniprot AS `{pathway}_uniprot` 
                ON (`{pathway}_targets_to_uniprot`.uniprot_id=`{pathway}_uniprot`.uniprot_id)
            INNER JOIN uniprot_to_pathway AS `{pathway}_uniprot_to_pathway` 
                ON (`{pathway}_targets_to_uniprot`.uniprot_id=`{pathway}_uniprot_to_pathway`.uniprot_id)
            INNER JOIN (
                SELECT pathway_id, COUNT(uniprot_id) AS `counts`
                FROM uniprot_to_pathway
                GROUP BY pathway_id
            ) AS `{pathway}_counts` ON `{pathway}_counts`.pathway_id=`{pathway}_uniprot_to_pathway`.pathway_id
            INNER JOIN pathway AS `{pathway}_pathway`
                ON (`{pathway}_uniprot_to_pathway`.pathway_id={pathway}_pathway.pathway_id)
        '''
        for pathway in pathway_names
    ))

    conditions = "\n".join((
        f'''
        AND `{pathway_name}_activity`.above_{threshold}=(1) 
        AND `{pathway_name}_pathway`.pathway_name="{pathway}"
        AND `{pathway_name}_pathway`.organism="{organism}"
        '''
        for pathway_name, pathway, organism in zip(pathway_names[1:], pathways[1:], organisms[1:])
    ))

    group_by = ",".join((
        f'''
            `{pathway_name}_pathway`.pathway_name,
            `{pathway_name}_pathway`.organism,
            `{pathway_name}_pathway`.pathway_url
       '''
        for pathway_name in pathway_names
    ))

    query = f'''
        SELECT c.compound_id, c.coconut_id AS `ID`, c.image_path AS `Image`, c.name AS `Molecule Name`, 
            c.formula AS `Molecular Formula`, c.clean_smiles AS `SMILES`,
        {columns}
        FROM compounds AS c
        {tables}
        WHERE `{pathway_names[0]}_activity`.above_{threshold}=(1)
        AND `{pathway_names[0]}_pathway`.pathway_name="{pathways[0]}"
        AND `{pathway_names[0]}_pathway`.organism="{organisms[0]}"
        {conditions}
        GROUP BY c.compound_id, c.coconut_id, c.name, c.formula, c.clean_smiles, {group_by}
        {f"LIMIT {limit}" if limit is not None else ""}
    '''

    records, cols = mysql_query(query, return_cols=True,)

    records = [
        record[1:] for record in records
    ]
    return records, cols[1:]

def query_reaction_hits(
    reactions,
    organisms,
    threshold=0,
    filter_pa_pi=True, 
    limit=None,
    existing_conn=None):
    assert filter_pa_pi

    if not isinstance(reactions, list):
        assert isinstance(reactions, str)
        reactions = [reactions]
    if not isinstance(organisms, list):
        assert isinstance(organisms, str)
        organisms = [organisms]

    reaction_names = sanitise_names(
        [   f"{reaction}_{organism}"
            for reaction, organism in zip(reactions, organisms)]
    )
    columns = ", ".join((
        f'''
        GROUP_CONCAT(DISTINCT(`{reaction}_target`.target_name) SEPARATOR '-') AS `{reaction} Target Names`, 
        COUNT(DISTINCT(`{reaction}_target`.target_name)) AS `{reaction} Number Targets`, 
        GROUP_CONCAT(DISTINCT(`{reaction}_uniprot`.acc) SEPARATOR '-') AS `{reaction} Uniprot ACCs`, 
        COUNT(DISTINCT(`{reaction}_uniprot`.acc)) AS `{reaction} Number Uniprot ACCs`,
        `{reaction}_counts`.counts AS `{reaction} Total Uniprot ACCs`,
        CAST(COUNT(DISTINCT(`{reaction}_uniprot`.acc)) / `{reaction}_counts`.counts AS CHAR)
            AS `{reaction} Uniprot Coverage`,    
        `{reaction}_reaction`.reaction_name AS `{reaction} Name`,
        `{reaction}_reaction`.organism AS `{reaction} Organism`,
        `{reaction}_reaction`.reaction_url AS `{reaction} URL`
        '''
        for reaction in reaction_names
    ))

    tables = "\n".join((
        f'''
            INNER JOIN activities AS `{reaction}_activity` 
                ON (c.compound_id=`{reaction}_activity`.compound_id)
            INNER JOIN targets AS `{reaction}_target` 
                ON (`{reaction}_activity`.target_id=`{reaction}_target`.target_id)
            INNER JOIN targets_to_uniprot AS `{reaction}_targets_to_uniprot` 
                ON (`{reaction}_activity`.target_id=`{reaction}_targets_to_uniprot`.target_id)
            INNER JOIN uniprot AS `{reaction}_uniprot` 
                ON (`{reaction}_targets_to_uniprot`.uniprot_id=`{reaction}_uniprot`.uniprot_id)
            INNER JOIN uniprot_to_reaction AS `{reaction}_uniprot_to_reaction` 
                ON (`{reaction}_targets_to_uniprot`.uniprot_id=`{reaction}_uniprot_to_reaction`.uniprot_id)
            INNER JOIN (
                SELECT reaction_id, COUNT(uniprot_id) AS `counts`
                FROM uniprot_to_reaction
                GROUP BY reaction_id
            ) AS `{reaction}_counts` ON `{reaction}_counts`.reaction_id=`{reaction}_uniprot_to_reaction`.reaction_id
            INNER JOIN reaction AS `{reaction}_reaction`
                ON (`{reaction}_uniprot_to_reaction`.reaction_id={reaction}_reaction.reaction_id)
        '''
        for reaction in reaction_names
    ))

    conditions = "\n".join((
        f'''
        AND `{reaction_name}_activity`.above_{threshold}=(1) 
        AND `{reaction_name}_reaction`.reaction_name="{reaction}"
        AND `{reaction_name}_reaction`.organism="{organism}"
        '''
        for reaction_name, reaction, organism in zip(reaction_names[1:], reactions[1:], organisms[1:])
    ))

    group_by = ",".join((
        f'''
            `{reaction_name}_reaction`.reaction_name,
            `{reaction_name}_reaction`.organism,
            `{reaction_name}_reaction`.reaction_url
       '''
        for reaction_name in reaction_names
    ))

    query = f'''
        SELECT c.compound_id, c.coconut_id AS `ID`, c.image_path AS `Image`, c.name AS `Molecule Name`, 
            c.formula AS `Molecular Formula`, c.clean_smiles AS `SMILES`,
        {columns}
        FROM compounds AS c
        {tables}
        WHERE `{reaction_names[0]}_activity`.above_{threshold}=(1)
        AND `{reaction_names[0]}_reaction`.reaction_name="{reactions[0]}"
        AND `{reaction_names[0]}_reaction`.organism="{organisms[0]}"
        {conditions}
        GROUP BY c.compound_id, c.coconut_id, c.name, c.formula, c.clean_smiles, {group_by}
        {f"LIMIT {limit}" if limit is not None else ""}
    '''

    records, cols = mysql_query(query, return_cols=True, )

    records = [
        record[1:] for record in records
    ]

    return records, cols[1:]

def get_multiple_compound_info(
    compounds=None, 
    columns=("coconut_id", "name", "formula", "clean_smiles"),
    limit=None):

    if compounds is not None:
        if not isinstance(compounds, str):
            if isinstance(compounds, list) or isinstance(compounds, set):
                compounds = tuple(compounds)
            assert isinstance(compounds, tuple)
            if len(compounds) == 1:
                compounds = compounds[0]

    compound_query = f'''
        SELECT {(", ".join(columns))}
        FROM compounds
        {f"WHERE coconut_id IN {compounds}" if isinstance(compounds, tuple)
        else f'WHERE coconut_id="{compounds}"' if isinstance(compounds, str) else ""}
        {f"LIMIT {limit}" if limit is not None else ""}
    '''

    return mysql_query(compound_query)

def get_coconut_compound_info_from_mongo( # use mongo
    compound_id, 
    projection={"_id": 0},
    compound_info_collection="uniqueNaturalProduct",
    filter_pa_pi=True,
    ):

    print ("querying COCONUT database for info",
        "about compound with compound id", compound_id)

    db = connect_to_mongodb()

    coconut_collection = db[compound_info_collection]

    compound_info = coconut_collection.find_one(
        {"coconut_id": compound_id},
        projection=projection) # get all info from mongo


    return compound_info

def get_all_activities_for_compound(
    coconut_id, 
    threshold=0,
    filter_pa_pi=True):
    print ("getting all activities for compound", coconut_id)

    categories = get_categories()
    category_targets = {category: get_targets_for_category(category)
        for category in categories}

    all_targets_query = f'''
        SELECT t.target_id, t.target_name, a.Pa, a.Pi, a.Pa-a.Pi
        FROM targets AS t
        INNER JOIN activities AS a ON (a.target_id=t.target_id)
        INNER JOIN compounds AS c ON (a.compound_id=c.compound_id)
        WHERE c.coconut_id='{coconut_id}'
        {f"AND a.above_{threshold}=(1)" if threshold > 0 and filter_pa_pi else ""}
        '''

    compounds_hits = mysql_query(all_targets_query)
    compound_hits = {i: (name, pa, pi, pa_pi) 
        for i, name, pa, pi, pa_pi in compounds_hits} #  convert to dict to avoid multiple queries

    category_activities = [("All Targets",
        [(name, pa, pi, pa_pi)
            for name, pa, pi, pa_pi in compound_hits.values()]
    )]

    category_activities += [(category.capitalize(),
        [(compound_hits[target_id][0], compound_hits[target_id][1], compound_hits[target_id][2], # name, pa, pi
            compound_hits[target_id][3]) # pa-pi
                for target_id in targets.values()
                    if target_id in compound_hits])
        for category, targets in category_targets.items()
    ]
    
    return category_activities

def draw_molecule(
    compound_id,
    smiles, 
    static_dir="static",
    output_dir="compound_images",
    ):
    output_dir = os.path.join(output_dir, f"{compound_id//1024}")
    os.makedirs(os.path.join(static_dir, output_dir), exist_ok=True)

    img_filename = os.path.join(output_dir, f"{compound_id}.png")
    img_full_path = os.path.join(static_dir, img_filename)
    if os.path.exists(img_full_path):
        print (img_full_path, "already exists")
        return img_filename
   
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        img = MolToImage(mol)
        img.save(img_full_path)
        print ("drawimg molecule to", img_full_path)
        return img_filename
    else:
        return None

def write_records_to_file(
    user_id,
    targets, 
    thresholds,
    records,
    root_dir="user_files",
    ):
    assert isinstance(records, pd.DataFrame)
    if "Image" in records.columns:
        del records["Image"]

    output_dir = os.path.join(root_dir, f"user={user_id}", "hit_records")
    os.makedirs(output_dir, exist_ok=True)

    targets = ",".join(map(lambda s: s.replace(" ", "_"), targets))
    thresholds= ",".join(map(str, thresholds))

    records_filename = os.path.join(output_dir,
        f'''targets={targets}-thresholds={thresholds}-results.csv''')
    print ("writing records to", records_filename)
    records.to_csv(records_filename)

    return records_filename

def write_smiles_to_file(
    user_id,
    targets,
    thresholds,
    smiles,
    root_dir="user_files",
    ):  

    output_dir = os.path.join(root_dir,f"user_id={user_id}", "smiles")
    os.makedirs(output_dir, exist_ok=True)

    targets = ",".join(map(lambda s: s.replace(" ", "_"), targets))
    thresholds= ",".join(map(str, thresholds))

    smiles_filename = os.path.join(output_dir,
        f'''targets={targets}-thresholds={thresholds}-hits.smi''')
    print ("writing smiles to", smiles_filename)

    write_smiles(smiles, smiles_filename)

    return smiles_filename
    
if __name__ == "__main__":

    from timeit import default_timer

    # "Stridor": 6773,
    # "Keratitis": 4205,
    # "Acidosis": 504,

    # start_time = default_timer()
    # hits, cols = query_target_hits([
    #     "Diarrhea", 
    #     "Yawning", 
    #     "Nausea",
    #     "Paralysis"
    # ], 900, 
    # filter_pa_pi=True,
    # # limit=100
    # )
    # print (len (hits))
    # print (default_timer() - start_time)
    # print (cols)
    # for hit in hits[:10]:
    #     print (hit)

    # compound_info, activities = get_compound_info("CNP0000002", filter_pa_pi=True)
    # all_activities = {c: a for c, a in get_all_activities_for_compound("CNP0000002")}

    # for row in all_activities["TOXICITY"]:
    #     print (row)

    # records = get_multiple_compound_info(columns=["coconut_id", "smiles"])
    # for record in records[:5]:
    #     print (record)

    # pathway_hits, cols = query_pathway_hits([
    #     "Zinc transporters", 
    #     "Xenobiotics",
    #     # "Disease",
    #     "Signaling by Erbb2",
    #     # "PI3K phosphorylates PIP2 to PIP3", # reaction
    #     # "Recruitment of PLCgamma to membrane"
    #     ],
    #     organism="Homo sapiens", 
    #     filter_pa_pi=True, threshold=950, 
    #     # limit=10
    # )

    # print (len(pathway_hits))
    # # for hit in pathway_hits[:1]:
    #     # print (hit)
    # print (cols)

    query = '''
        SELECT compound_id, clean_smiles
        FROM compounds
    '''
    records = mysql_query(query)

    for _id, smiles in records:
        draw_molecule(_id, smiles)