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
from utils.mysql_utils import mysql_query, mysql_create_table, mysql_insert_many, connect_to_mysqldb
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

    target_names = [re.sub(r"( |-|\(|\))", "_", target)
        for target in targets]
    columns = ", ".join((
        f'''
        `{target}_activity`.Pa AS `{target}-Pa`, `{target}_activity`.Pi AS `{target}-Pi`, 
            `{target}_activity`.Pa-`{target}_activity`.Pi AS `{target}-Pa-Pi`
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
        AND `{target_name}_target`.target_name='{target}'
        '''
        for target_name, target, threshold in zip(target_names[1:], targets[1:], thresholds[1:])      
    ))
    query = f'''
        SELECT c.compound_id, c.coconut_id AS `ID`, c.name AS `Molecule Name`, 
            c.formula AS `Molecular Formula`, c.clean_smiles AS `SMILES`,
        {columns}
        FROM compounds AS c
        {tables}
        WHERE `{target_names[0]}_activity`.above_{thresholds[0]}=(1)
        AND `{target_names[0]}_target`.target_name='{targets[0]}'
        {conditions}
        {f"LIMIT {limit}" if limit is not None else ""}
    '''

    hits, cols = mysql_query(query, return_cols=True)

    records = [
        record[1:] for record in hits
    ]

    return records, cols[1:]

def query_pathway_hits(pathways,
    organism=None, threshold=0, filter_pa_pi=True, 
    include_targets=True,
    limit=None,
    existing_conn=None):
    assert filter_pa_pi
    # assert threshold in (900, )

    # if existing_conn is None:
    #     existing_conn = connect_to_mysqldb()

    if not isinstance(pathways, list):
        assert isinstance(pathways, str)
        pathways = [pathways]

    pathway_names = [re.sub(r"( |-|\(|\)|/)", "_", pathway)
        for pathway in pathways]

    columns = ", ".join((
        f'''
        GROUP_CONCAT(DISTINCT(`{pathway}_target`.target_name)) AS `{pathway} Target Names`, 
        COUNT(DISTINCT(`{pathway}_target`.target_name)) AS `{pathway} Number Targets`, 
        GROUP_CONCAT(DISTINCT(`{pathway}_uniprot`.acc)) AS `{pathway} Uniprot ACCs`, 
        COUNT(DISTINCT(`{pathway}_uniprot`.acc)) AS `{pathway} Number Uniprot ACCs`, 
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
            INNER JOIN pathway AS `{pathway}_pathway`
                ON (`{pathway}_uniprot_to_pathway`.pathway_id={pathway}_pathway.pathway_id)
        '''
        for pathway in pathway_names
    ))

    conditions = "\n".join((
        f'''
        AND `{pathway_name}_activity`.above_{threshold}=(1) 
        AND `{pathway_name}_pathway`.pathway_name='{pathway}'
        AND `{pathway_name}_pathway`.organism='{organism}'
        '''
        for pathway_name, pathway in zip(pathway_names[1:], pathways[1:])
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
        SELECT c.compound_id, c.coconut_id AS `ID`, c.name AS `Molecule Name`, 
            c.formula AS `Molecular Formula`, c.clean_smiles AS `SMILES`,
        {columns}
        FROM compounds AS c
        {tables}
        WHERE `{pathway_names[0]}_activity`.above_{threshold}=(1)
        AND `{pathway_names[0]}_pathway`.pathway_name='{pathways[0]}'
        AND `{pathway_names[0]}_pathway`.organism='{organism}'
        {conditions}
        GROUP BY c.compound_id, c.coconut_id, c.name, c.formula, c.clean_smiles, {group_by}
        {f"LIMIT {limit}" if limit is not None else ""}
    '''

    # # get unique compounds that hit all pathways

    #  # initial query
    # pathway_name = pathway_names[0] 
    # query = f'''
    #     SELECT DISTINCT a.compound_id
    #     FROM activities AS a
    #     INNER JOIN targets_to_uniprot AS tu ON (a.target_id=tu.target_id)
    #     INNER JOIN uniprot_to_pathway AS up ON (tu.uniprot_id=up.uniprot_id)
    #     INNER JOIN pathway AS p ON (up.pathway_id=p.pathway_id)
    #     WHERE p.pathway_name="{pathway_name}"
    #     {f'AND p.organism="{organism}"' if organism is not None else ""}
    #     AND a.above_{threshold}=(1)
    #     {f"LIMIT {limit}" if limit is not None else ""}
    # '''
    # # save in case of need
    # # {f"AND a.Pa>{threshold}" if threshold>0 else ""}
    # # {"AND a.Pa>a.Pi" if filter_pa_pi else ""}

    # for pathway_name in pathway_names[1:]:

    #     query = f'''
    #         SELECT DISTINCT q.compound_id
    #         FROM ({query}) as q
    #         INNER JOIN activities AS a ON (q.compound_id=a.compound_id)
    #         INNER JOIN targets_to_uniprot AS tu ON (a.target_id=tu.target_id)
    #         INNER JOIN uniprot_to_pathway AS up ON (tu.uniprot_id=up.uniprot_id)
    #         INNER JOIN pathway AS p ON up.pathway_id=p.pathway_id
    #         WHERE p.pathway_name="{pathway_name}"
    #         {f'AND p.organism="{organism}"' if organism is not None else ""}
    #         AND a.above_{threshold}=(1)
    #         {f"LIMIT {limit}" if limit is not None else ""}
    #     '''

    # compound_hits = mysql_query(query, existing_conn=existing_conn)
    # compound_hits = [hit[0] for hit in compound_hits]
    # num_compound_hits = len(compound_hits)

    # print ("number of unique compounds that hit all pathways:",  num_compound_hits)
    # if num_compound_hits == 0:
    #     return [], None

    # query = f'''
    #     SELECT DISTINCT c.coconut_id AS 'ID', c.name AS 'Molecule Name', c.formula AS 'Molecular Formula',
    #         {f"t.target_name AS 'Target Name', u.acc AS 'Target UNIPROT ACC', a.Pa AS 'Pa', a.Pi AS 'Pi',"
    #         if include_targets else ""}
    #         up.evidence AS 'Evidence', p.pathway_name as 'Pathway Name', p.organism AS 'Organism', p.pathway_url AS 'Pathway URL'
    #     FROM compounds AS c
    #     INNER JOIN activities AS a ON (c.compound_id=a.compound_id)
    #     INNER JOIN targets_to_uniprot AS tu ON (a.target_id=tu.target_id)
    #     {"INNER JOIN targets AS t ON (a.target_id=t.target_id)" 
    #         if include_targets else ""}
    #     {"INNER JOIN uniprot AS u ON (tu.uniprot_id=u.uniprot_id)"
    #         if include_targets else ""}
    #     INNER JOIN uniprot_to_pathway AS up ON (tu.uniprot_id=up.uniprot_id)
    #     INNER JOIN pathway AS p ON (up.pathway_id=p.pathway_id)
    #     WHERE {f"p.pathway_name IN {tuple(pathway_names)}" if len(pathway_names)>1
    #         else f'p.pathway_name="{pathway_names[0]}"'}
    #     {f'AND p.organism="{organism}"' if organism is not None else ""}
    #     AND a.above_{threshold}=(1)
    #     AND c.compound_id IN {tuple(compound_hits)}
    # ''' 

    records, cols = mysql_query(query, return_cols=True, 
        # existing_conn=existing_conn
    )

    records = [
        record[1:] for record in records
    ]

    return records, cols[1:]

def query_reaction_hits(reactions,
    organism=None, threshold=0, filter_pa_pi=True, 
    include_targets=True,
    limit=None,
    existing_conn=None):
    assert filter_pa_pi
    # assert threshold in (900, )

    if not isinstance(reactions, list):
        assert isinstance(reactions, str)
        reactions = [reactions]

    reaction_names = [re.sub(r"( |-|\(|\)|/)", "_", reaction)
        for reaction in reactions]

    columns = ", ".join((
        f'''
        GROUP_CONCAT(DISTINCT(`{reaction}_target`.target_name)) AS `{reaction} Target Names`, 
        COUNT(DISTINCT(`{reaction}_target`.target_name)) AS `{reaction} Number Targets`, 
        GROUP_CONCAT(DISTINCT(`{reaction}_uniprot`.acc)) AS `{reaction} Uniprot ACCs`, 
        COUNT(DISTINCT(`{reaction}_uniprot`.acc)) AS `{reaction} Number Uniprot ACCs`, 
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
            INNER JOIN reaction AS `{reaction}_reaction`
                ON (`{reaction}_uniprot_to_reaction`.reaction_id={reaction}_reaction.reaction_id)
        '''
        for reaction in reaction_names
    ))

    conditions = "\n".join((
        f'''
        AND `{reaction_name}_activity`.above_{threshold}=(1) 
        AND `{reaction_name}_reaction`.reaction_name='{reaction}'
        AND `{reaction_name}_reaction`.organism='{organism}'
        '''
        for reaction_name, reaction in zip(reaction_names[1:], reactions[1:])
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
        SELECT c.compound_id, c.coconut_id AS `ID`, c.name AS `Molecule Name`, 
            c.formula AS `Molecular Formula`, c.clean_smiles AS `SMILES`,
        {columns}
        FROM compounds AS c
        {tables}
        WHERE `{reaction_names[0]}_activity`.above_{threshold}=(1)
        AND `{reaction_names[0]}_reaction`.reaction_name='{reactions[0]}'
        AND `{reaction_names[0]}_reaction`.organism='{organism}'
        {conditions}
        GROUP BY c.compound_id, c.coconut_id, c.name, c.formula, c.clean_smiles, {group_by}
        {f"LIMIT {limit}" if limit is not None else ""}
    '''


    # # get unique compounds that hit all reactions

    #  # initial query
    # reaction_name = reaction_names[0] 
   
    # query = f'''
    #     SELECT DISTINCT a.compound_id
    #     FROM activities AS a
    #     INNER JOIN targets_to_uniprot AS tu ON (a.target_id=tu.target_id)
    #     INNER JOIN uniprot_to_reaction AS ur ON (tu.uniprot_id=ur.uniprot_id)
    #     INNER JOIN reaction AS r ON (ur.reaction_id=r.reaction_id)
    #     WHERE r.reaction_name="{reaction_name}"
    #     {f'AND r.organism="{organism}"' if organism is not None else ""}
    #     AND a.above_{threshold}=(1)
    #     {f"LIMIT {limit}" if limit is not None else ""}
    # '''

    # for reaction_name in reaction_names[1:]:
    #     query = f'''
    #         SELECT DISTINCT a.compound_id
    #         FROM ({query}) as q
    #         INNER JOIN activities AS a ON (q.compound_id=a.compound_id)
    #         INNER JOIN targets_to_uniprot AS tu ON (a.target_id=tu.target_id)
    #         INNER JOIN uniprot_to_reaction AS ur ON (tu.uniprot_id=ur.uniprot_id)
    #         INNER JOIN reaction AS r ON ( ur.reaction_id=r.reaction_id)
    #         WHERE r.reaction_name="{reaction_name}"
    #         {f'AND r.organism="{organism}"' if organism is not None else ""}
    #         AND a.above_{threshold}=(1)
    #         {f"LIMIT {limit}" if limit is not None else ""}
    #     '''

    # compound_hits = mysql_query(query, existing_conn=existing_conn)
    # compound_hits = [hit[0] for hit in compound_hits]
    # num_compound_hits = len(compound_hits)

    # print ("number of unique compounds that hit all reactions:",  num_compound_hits)
    # if num_compound_hits == 0:
    #     return [], None

    # query = f'''
    #     SELECT DISTINCT c.coconut_id AS 'ID', c.name AS 'Molecule Name', c.formula AS 'Molecular Formula',
    #         {f"t.target_name AS 'Target Name', u.acc AS 'Target UNIPROT ACC', a.Pa AS 'Pa', a.Pi AS 'Pi',"
    #         if include_targets else ""}
    #         ur.evidence AS 'Evidence', r.reaction_name AS 'Reaction Name', r.organism AS 'Organism', r.reaction_url AS 'Reaction URL'
    #     FROM compounds AS c
    #     INNER JOIN activities AS a ON (c.compound_id=a.compound_id)
    #     INNER JOIN targets_to_uniprot AS tu ON (a.target_id=tu.target_id)
    #     {"INNER JOIN targets AS t ON (a.target_id=t.target_id)" 
    #         if include_targets else ""}
    #     {"INNER JOIN uniprot AS u ON (tu.uniprot_id=u.uniprot_id)"
    #         if include_targets else ""}
    #     INNER JOIN uniprot_to_reaction AS ur ON (tu.uniprot_id=ur.uniprot_id)
    #     INNER JOIN reaction AS r ON (ur.reaction_id=r.reaction_id)
    #     WHERE {f"r.reaction_name in {tuple(reaction_names)}" if len(reaction_names)>1
    #         else f'r.reaction_name="{reaction_names[0]}"'}
    #     {f'AND r.organism="{organism}"' if organism is not None else ""}
    #     AND a.above_{threshold}=(1)
    #     AND c.compound_id IN {tuple(compound_hits)}
    # ''' 

    records, cols = mysql_query(query, return_cols=True, 
        # existing_conn=existing_conn
    )

    records = [
        record[1:] for record in records
    ]

    return records, cols[1:]


def get_multiple_compound_info(
    compounds=None, 
    columns=("coconut_id", "name", "formula", "clean_smiles")):

    if compounds is not None:
        if not isinstance(compounds, tuple):
            compounds = tuple(compounds)

    compound_query = f'''
        SELECT {(", ".join(columns))}
        FROM compounds
        {f"WHERE coconut_id IN {compounds}" if compounds is not None else ""}
    '''

    records = mysql_query(compound_query)

    return records

def get_compound_info( # used mongo
    compound_id, 
    projection={"_id": 0},
    get_activities=True,
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

    if not get_activities:
        return compound_info

    pass_activities = get_all_activities_for_compound(
        compound_id, 
        filter_pa_pi=filter_pa_pi)

    return compound_info, pass_activities


def get_all_activities_for_compound(
    coconut_id, 
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
        {"AND a.Pa>a.Pi" if filter_pa_pi else ""}
        '''

    compounds_hits = mysql_query(all_targets_query)
    compound_hits = {i: (name, pa, pi, pa_pi) 
        for i, name, pa, pi, pa_pi in compounds_hits} #  convert to dict to avoid multiple queries

    category_activities = [("ALL",
        [(name, pa, pi, pa_pi)
            for name, pa, pi, pa_pi in compound_hits.values()]
    )]

    category_activities += [(category,
        [(compound_hits[target_id][0], compound_hits[target_id][1], compound_hits[target_id][2], # name, pa, pi
            compound_hits[target_id][3]) # pa-pi
                for target_id in targets.values()
                    if target_id in compound_hits])
        for category, targets in category_targets.items()
    ]
    
    return category_activities
    

def draw_molecule(smiles, 
    static_dir="natural_products/static",
    img_filename="natural_products/temp.png", 
    ):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        img = MolToImage(mol)

        img.save(os.path.join(static_dir,
            img_filename))

        return img_filename
    else:
        return None

def write_records_to_file(
    user_id,
    targets, 
    thresholds,
    records,
    static_dir="natural_products/static/natural_products",
    output_dir="results",
    ):
    assert isinstance(records, pd.DataFrame)

    output_dir = os.path.join(static_dir, output_dir, f"user={user_id}")
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
    static_dir="natural_products/static",
    output_dir="smiles",
    ):  

    output_dir = os.path.join(static_dir, output_dir, f"user_id={user_id}")
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

    pathway_hits, cols = query_pathway_hits([
        "Zinc transporters", 
        "Xenobiotics",
        # "Disease",
        "Signaling by Erbb2",
        # "PI3K phosphorylates PIP2 to PIP3", # reaction
        # "Recruitment of PLCgamma to membrane"
        ],
        organism="Homo sapiens", 
        filter_pa_pi=True, threshold=950, 
        # limit=10
    )

    print (len(pathway_hits))
    # for hit in pathway_hits[:1]:
        # print (hit)
    print (cols)