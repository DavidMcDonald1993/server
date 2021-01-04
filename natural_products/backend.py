import os
import sys
import os.path
sys.path.insert(1, 
    os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))

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

def query_target_hits(
    targets, 
    thresholds,
    filter_pa_pi=True,
    ):

    assert isinstance(targets, list)
    if not isinstance(thresholds, list):
        assert isinstance(thresholds, int), thresholds
        thresholds = [thresholds]
    num_targets = len(targets)
    if len(thresholds) < num_targets:
        thresholds = thresholds * num_targets

    print ("querying PASS database for compounds",
        "that hit targets", targets,
         "with Pa greater than or equal to",
        "thresholds", thresholds)

    # initial query
    target = targets[0] 
    threshold = thresholds[0]
    query = f'''SELECT c.compound_id, c.coconut_id, c.name, c.formula,
        a.Pa AS '{target}-Pa', a.Pi AS '{target}-Pi', a.Pa-a.Pi AS '{target}-Pa-Pi'
        FROM compounds AS c, activities AS a, targets AS t
        WHERE t.target_name='{target}' AND t.target_id=a.target_id
        AND a.compound_id=c.compound_id
        AND a.Pa>{threshold}
        {"AND a.Pa>a.Pi" if filter_pa_pi else ""}
    '''

    for target, threshold in zip(targets[1:], thresholds[1:]):
        query = f'''
        SELECT q.*, a.Pa AS '{target}-Pa', a.Pi AS '{target}-Pi', a.Pa-a.Pi AS '{target}-Pa-Pi'
        FROM ({query}) AS q, activities AS a, targets AS t
        WHERE t.target_name='{target}' AND t.target_id=a.target_id
        AND a.compound_id=q.compound_id
        AND a.Pa>{threshold}
        {"AND a.Pa>a.Pi" if filter_pa_pi else ""}
        '''

    hits =  mysql_query(query)

    # n_hits = len(hits)
    # print ("built records, number of hits", n_hits)

    records = [
        record[1:] for record in hits
    ]

    return records

def query_pathway_hits(pathway_names,
    organism=None, threshold=0, filter_pa_pi=True, 
    include_targets=True,
    limit=None,
    existing_conn=None):

    if existing_conn is None:
        existing_conn = connect_to_mysqldb()

    if not isinstance(pathway_names, list):
        assert isinstance(pathway_names, str)
        pathway_names = [pathway_names]

     # initial query
    pathway_name = pathway_names[0] 

    # query = f'''
    #     SELECT DISTINCT c.compound_id, c.coconut_id, c.name, c.formula,
    #         {f"t.target_name AS '{pathway_name}-target_name', u.acc AS '{pathway_name}-acc', a.Pa AS '{pathway_name}-Pa', a.Pi AS '{pathway_name}-Pi',"
    #         if include_targets else ""}
    #         p.pathway_name AS '{pathway_name}-pathway_name', p.organism AS '{pathway_name}-organism'
    #     FROM compounds AS c, activities AS a,
    #         targets_to_uniprot AS tu,
    #         {"targets AS t, uniprot AS u," 
    #             if include_targets else ""}
    #         uniprot_to_pathway AS up,
    #         pathway AS p
    #     WHERE c.compound_id=a.compound_id
    #     AND a.target_id=tu.target_id
    #     {"AND a.target_id=t.target_id AND tu.uniprot_id=u.uniprot_id" 
    #         if include_targets else ""}
    #     AND tu.uniprot_id=up.uniprot_id
    #     AND up.pathway_id=p.pathway_id
    #     AND p.pathway_name="{pathway_name}"
    #     {f"AND a.Pa>{threshold}" if threshold>0 else ""}
    #     {"AND a.Pa>a.Pi" if filter_pa_pi else ""}
    #     {f'AND p.organism="{organism}"' if organism is not None else ""}
    #     {f"LIMIT {limit}" if limit is not None else ""}
    # '''

    # get unique compounds that hit all pathways
    query = f'''
        SELECT DISTINCT a.compound_id
        FROM activities AS a,
            targets_to_uniprot AS tu,
            uniprot_to_pathway AS up,
            pathway AS p
        WHERE a.target_id=tu.target_id
        AND tu.uniprot_id=up.uniprot_id
        AND up.pathway_id=p.pathway_id
        AND p.pathway_name="{pathway_name}"
        {f'AND p.organism="{organism}"' if organism is not None else ""}
        {f"AND a.Pa>{threshold}" if threshold>0 else ""}
        {"AND a.Pa>a.Pi" if filter_pa_pi else ""}
        {f"LIMIT {limit}" if limit is not None else ""}
    '''

    for pathway_name in pathway_names[1:]:
        # query = f'''
        #     SELECT DISTINCT q.*,
        #         {f"t.target_name AS '{pathway_name}-target_name',  u.acc AS '{pathway_name}-acc', a.Pa AS '{pathway_name}-Pa', a.Pi AS '{pathway_name}-Pi',"
        #         if include_targets else ""}
        #         p.pathway_name AS '{pathway_name}-pathway_name', p.organism AS '{pathway_name}-organism'
        #     FROM ({query}) AS q, activities AS a,
        #         targets_to_uniprot AS tu,
        #         {"targets AS t, uniprot AS u," 
        #             if include_targets else ""}
        #         uniprot_to_pathway AS up,
        #         pathway AS p
        #     WHERE q.compound_id=a.compound_id
        #     AND a.target_id=tu.target_id
        #     {"AND a.target_id=t.target_id AND tu.uniprot_id=u.uniprot_id" 
        #         if include_targets else ""}
        #     AND tu.uniprot_id=up.uniprot_id
        #     AND up.pathway_id=p.pathway_id
        #     AND p.pathway_name="{pathway_name}"
        #     {f"AND a.Pa>{threshold}" if threshold>0 else ""}
        #     {"AND a.Pa>a.Pi" if filter_pa_pi else ""}
        #     {f'AND p.organism="{organism}"' if organism is not None else ""}
        #     {f"LIMIT {limit}" if limit is not None else ""}
        # '''

        query = f'''
            SELECT DISTINCT a.compound_id
            FROM ({query}) as q, activities AS a,
                targets_to_uniprot AS tu,
                uniprot_to_pathway AS up,
                pathway AS p
            WHERE q.compound_id=a.compound_id
            AND a.target_id=tu.target_id
            AND tu.uniprot_id=up.uniprot_id
            AND up.pathway_id=p.pathway_id
            AND p.pathway_name="{pathway_name}"
            {f'AND p.organism="{organism}"' if organism is not None else ""}
            {f"AND a.Pa>{threshold}" if threshold>0 else ""}
            {"AND a.Pa>a.Pi" if filter_pa_pi else ""}
            {f"LIMIT {limit}" if limit is not None else ""}
        '''

    hits = mysql_query(query, existing_conn=existing_conn)
    hits = [hit[0] for hit in hits]

    print ("number of unique compounds that hit all pathways:", len(hits))

    # records = []
    # for pathway_name in pathway_names:

    query = f'''
        SELECT DISTINCT c.coconut_id, c.name, c.formula,
            {f"t.target_name, u.acc, a.Pa, a.Pi,"
            if include_targets else ""}
            p.pathway_name, p.organism
        FROM compounds AS c, activities AS a,
            targets_to_uniprot AS tu,
            {"targets AS t, uniprot AS u," 
                if include_targets else ""}
            uniprot_to_pathway AS up,
            pathway AS p
        WHERE c.compound_id=a.compound_id
        AND a.target_id=tu.target_id
        {"AND a.target_id=t.target_id AND tu.uniprot_id=u.uniprot_id" 
            if include_targets else ""}
        AND tu.uniprot_id=up.uniprot_id
        AND up.pathway_id=p.pathway_id
        AND p.pathway_name in {tuple(pathway_names)}
        {f'AND p.organism="{organism}"' if organism is not None else ""}
        {f"AND a.Pa>{threshold}" if threshold>0 else ""}
        {"AND a.Pa>a.Pi" if filter_pa_pi else ""}
        AND c.compound_id IN {tuple(hits)}
    ''' 

        # records.append((pathway_name, mysql_query(query, existing_conn=existing_conn)))

    records = mysql_query(query, existing_conn=existing_conn)

    # return [hit[1:] for hit in hits]
    return records

def query_reaction_hits(reaction_names,
    organism=None, threshold=0, filter_pa_pi=True, 
    include_targets=True,
    limit=None,
    existing_conn=None):

    if not isinstance(reaction_names, list):
        assert isinstance(reaction_names, str)
        reaction_names = [reaction_names]

     # initial query
    reaction_name = reaction_names[0] 

    # query = f'''
    #     SELECT DISTINCT c.coconut_id, c.name, c.formula
    #         {", t.target_name AS {reaction_name}-target_name, u.acc" if include_targets else ""}
    #     FROM compounds AS c, activities AS a,
    #         targets_to_uniprot as tu,
    #         {"targets as t, uniprot AS u," 
    #             if include_targets else ""}
    #         uniprot_to_reaction AS ur,
    #         reaction AS r
    #     WHERE c.compound_id=a.compound_id
    #     AND a.target_id=tu.target_id
    #     {"AND a.target_id=t.target_id AND tu.uniprot_id=u.uniprot_id" 
    #         if include_targets else ""}
    #     AND tu.uniprot_id=ur.uniprot_id
    #     AND ur.reaction_id=r.reaction_id
    #     AND r.reaction_name="{reaction_name}"
    # '''
    # if threshold>0:
    #     query += f" AND a.Pa>{threshold}"
    # if filter_pa_pi:
    #     query += " AND a.Pa>a.Pi"
    # if organism is not None:
    #     query += f" AND r.organism=\"{organism}\""
   
    # get unique compounds that hit all reactions
    query = f'''
        SELECT DISTINCT a.compound_id
        FROM activities AS a,
            targets_to_uniprot AS tu,
            uniprot_to_reaction AS ur,
            reaction AS r
        WHERE a.target_id=tu.target_id
        AND tu.uniprot_id=ur.uniprot_id
        AND ur.reaction_id=r.reaction_id
        AND r.reaction_name="{reaction_name}"
        {f'AND r.organism="{organism}"' if organism is not None else ""}
        {f"AND a.Pa>{threshold}" if threshold>0 else ""}
        {"AND a.Pa>a.Pi" if filter_pa_pi else ""}
        {f"LIMIT {limit}" if limit is not None else ""}
    '''

    for reaction_name in reaction_names[1:]:
        query = f'''
            SELECT DISTINCT a.compound_id
            FROM ({query}) as q, activities AS a,
                targets_to_uniprot AS tu,
                uniprot_to_reaction AS ur,
                reaction AS r
            WHERE q.compound_id=a.compound_id
            AND a.target_id=tu.target_id
            AND tu.uniprot_id=ur.uniprot_id
            AND ur.reaction_id=r.reaction_id
            AND r.reaction_name="{reation_name}"
            {f'AND r.organism="{organism}"' if organism is not None else ""}
            {f"AND a.Pa>{threshold}" if threshold>0 else ""}
            {"AND a.Pa>a.Pi" if filter_pa_pi else ""}
            {f"LIMIT {limit}" if limit is not None else ""}
        '''

    hits = mysql_query(query, existing_conn=existing_conn)
    hits = [hit[0] for hit in hits]

    print ("number of unique compounds that hit all reactions:", len(hits))

    query = f'''
        SELECT DISTINCT c.coconut_id, c.name, c.formula,
            {f"t.target_name, u.acc, a.Pa, a.Pi,"
            if include_targets else ""}
            r.reaction_name, r.organism
        FROM compounds AS c, activities AS a,
            targets_to_uniprot AS tu,
            {"targets AS t, uniprot AS u," 
                if include_targets else ""}
            uniprot_to_reaction AS ur,
            reaction AS r
        WHERE c.compound_id=a.compound_id
        AND a.target_id=tu.target_id
        {"AND a.target_id=t.target_id AND tu.uniprot_id=u.uniprot_id" 
            if include_targets else ""}
        AND tu.uniprot_id=ur.uniprot_id
        AND ur.reaction_id=r.reaction_id
        AND r.reaction_name in {tuple(reaction_names)}
        {f'AND r.organism="{organism}"' if organism is not None else ""}
        {f"AND a.Pa>{threshold}" if threshold>0 else ""}
        {"AND a.Pa>a.Pi" if filter_pa_pi else ""}
        AND c.compound_id IN {tuple(hits)}
    ''' 

    records = mysql_query(query, existing_conn=existing_conn)

    return records


def get_multiple_compound_info(
    compounds=None, 
    columns=("coconut_id", "name", "formula", "clean_smiles")):

    if compounds is not None:
        if not isinstance(compounds, tuple):
            compounds = tuple(compounds)

    compound_query = f'''
        SELECT {(", ".join(columns))}
        FROM compounds
        {"WHERE coconut_id IN {compounds}" if compounds is not None else ""}
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
        FROM targets AS t, activities AS a, compounds AS c
        WHERE c.coconut_id='{coconut_id}'
        AND a.compound_id=c.compound_id AND a.target_id=t.target_id
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
    # mol = standardise_smi(smiles)
    if mol is not None:
        img = MolToImage(mol)

        img.save(os.path.join(static_dir,
            img_filename))

        return img_filename
    else:
        return None

def write_smiles_to_file(
    username,
    smiles,
    static_dir="natural_products/static",
    output_dir="smiles",
    ):  

    output_dir = os.path.join(static_dir, output_dir, username)

    os.makedirs(output_dir, exist_ok=True)

    smiles_filename = os.path.join(output_dir,
        "results.smi")
    print ("writing smiles to", smiles_filename)

    # columns are ["coconut_id",  "molecule_name", "molecular_formula", "smiles", "Pa", "Pi", "Pa-Pi", ...]
    # smiles = [(record[0], record[3]) for record in records]

    write_smiles(smiles, smiles_filename)

    return smiles_filename

def write_records_to_file(
    username,
    targets, 
    thresholds,
    records,
    static_dir="natural_products/static",
    output_dir="results",
    ):

    output_dir = os.path.join(static_dir, output_dir, username)

    os.makedirs(output_dir, exist_ok=True)

    columns = ["coconut_id",  "molecule_name", "molecular_formula", ] #"smiles", ]
    for target in targets:
        columns += [f"{target}-Pa", f"{target}-Pi", f"{target}-Pa-Pi" ]

    records = pd.DataFrame.from_records(records, 
        columns=columns)

    records_filename = os.path.join(output_dir,
        f"{targets}-{thresholds}-results.csv")
    print ("writing records to", records_filename)
    records.to_csv(records_filename)

    return records_filename
    
if __name__ == "__main__":

    # hits = query_target_hits(["Diarrhea", "Yawning", "Nausea"], 500, filter_pa_pi=False)
    # print (hits)

    # compound_info, activities = get_compound_info("CNP0000002", filter_pa_pi=True)
    # all_activities = {c: a for c, a in get_all_activities_for_compound("CNP0000002")}

    # for row in all_activities["TOXICITY"]:
    #     print (row)

    # records = get_multiple_compound_info(columns=["coconut_id", "smiles"])
    # for record in records[:5]:
    #     print (record)

    pathway_hits = query_pathway_hits(pathway_names=[
        "Zinc transporters", 
        "Xenobiotics",
        "Disease",
        ],
        organism="Homo sapiens", filter_pa_pi=True, threshold=950, limit=10000)

    print (len(pathway_hits))
    # for hit in pathway_hits[:10]:
        # print (hit)
    # for pathway, hits in pathway_hits:
        # print (pathway)
        # for record in hits[:10]:
            # print (record)
        # print ()