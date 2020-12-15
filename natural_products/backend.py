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
from utils.pass_utils import get_categories, get_all_targets, get_targets_for_category

# PASS_COLLECTION = "PASS_ALL"

# target_categories = [
#     "EFFECTS",
#     "MECHANISMS",
#     'TOXICITY', 
#     'ANTITARGETS',
#     'METABOLISM', 
#     'GENE_EXPRESSION', 
#     'TRANSPORTERS'
# ]

# def get_categories():
#     return target_categories

# def get_targets_for_category(category,
#     collection=PASS_COLLECTION):
#     assert category in target_categories

#     print ("identifying targets for category", category)

#     db = connect_to_mongodb()
#     pass_collection = db[collection]

#     query = {"coconut_id": "CNP0000002"}
#     projection = {"_id": 0, "PASS_" + category: 1}

#     record = pass_collection.find_one(query, projection)

#     return list(record["PASS_" + category].keys())

def query_target_hits(
    # category,
    target, 
    threshold,
    filter_pa_pi=True,
    # collection=PASS_COLLECTION
    ):
    # assert category in target_categories
    # category = "PASS_" + category

    print ("querying PASS database for compounds",
        "that hit target", target, #"in category", category,
         "with Pa greater than or equal to",
        "threshold", threshold)

    target_ids = get_all_targets()

    target_id = target_ids[target]

    query = f"SELECT coconut_id, Pa, Pi from `{target_id}` where Pa>{threshold}"
    if filter_pa_pi:
        query += " AND Pa>Pi"
    
    print ("using query string", query)

    hits = mysql_query(query)

    # db = connect_to_mongodb()
    # pass_collection = db[collection]

    # category_target = "{}.{}".format(category, target)

    # query = [{"{}.Pa".format(category_target): {"$gte": threshold}, }]
    # if filter_pa_pi:
    #     query.append( 
    #         {"$expr": { "$gt": 
    #             ["${}.Pa".format(category_target), "${}.Pi".format(category_target)] } } )
    
    # projection = {
    #     "_id": 0, 
    #     "coconut_id": 1, 
    #     # "name": 1,
    #     # "molecular_formula": 1, 
    #     # "SMILES": 1, # get these from compound database
    #     category_target: 1
    # }

    # query = {"$and": query}

    # print ("filtering with", query)
    # print ("showing", projection)

    # cursor = pass_collection.find(query, projection).sort("coconut_id", 1).limit(100)
    
    # print ("iterating over records")

    # records = [
    #     (record["coconut_id"], 
    #         # (urlparse.unquote(record["name"]).capitalize() 
    #             # if record["name"] is not None else "<NO NAME>"),
    #         # record["molecular_formula"], 
    #         # record["SMILES"], 
    #         record[category][target]["Pa"], record[category][target]["Pi"], 
    #         record[category][target]["Pa"] - record[category][target]["Pi"])
    #     for record in cursor
    # ]

    # db.client.close()

    n_hits = len(hits)

    print ("built records, number of hits", n_hits)

    # get additional compound information from compound databaase

    print ("getting info about", n_hits, "compounds")
    coconut_ids = [hit[0] for hit in hits]
    compound_info = get_multiple_compound_info(
        coconut_ids)

    # compound_info_records = [
    #     get_compound_info(coconut_id, 
    #         projection={
    #             "_id": 0, 
    #             "coconut_id": 1, 
    #             "name": 1,
    #             "molecular_formula": 1, 
    #             "clean_smiles": 1, 
    #         },
    #         get_activities=False)

    #     for coconut_id in coconut_ids
    # ]

    # assert len(compound_info_records) == n_compounds

    # for record, compound_info in zip(records, compound_info_records):
        # assert record[0] == compound_info[0]

    records = [
            (hit[0], # coconut_id
                compound_info[hit[0]][0], # name
                compound_info[hit[0]][1], # molecular_formula 
                compound_info[hit[0]][2], # clean_smiles
                hit[1], # Pa
                hit[2], # Pi
                hit[1] - hit[2]) # Pa - Pi
        for hit in hits
        if hit[0] in compound_info
    ]

    return records[:10]


def get_multiple_compound_info(compounds=None,    
    compound_info_collection="uniqueNaturalProduct",):

    print ("querying COCONUT database for info",
        "about all compounds")

    db = connect_to_mongodb()

    coconut_collection = db[compound_info_collection]

    if compounds is None: # get all
        query = {}
    else:
        if not isinstance(compounds, list):
            compounds = list(compounds)
        query = {"coconut_id": {"$in": compounds}}

    projection = {
        "_id": 0, 
        "coconut_id": 1, 
        "name": 1,
        "molecular_formula": 1, 
        "clean_smiles": 1
    }

    cursor = coconut_collection.find(query, projection).sort("coconut_id", 1)
    # cursor.batch_size(1000000)

    print ("iterating over query")

    records = {
        record["coconut_id"]:
            ((urlparse.unquote(record["name"]).capitalize() 
                    if "name" in record else "<NO NAME>"),
            record["molecular_formula"],
            record["clean_smiles"])
        for record in cursor
    }

    db.client.close()

    print ("built records", len(records))

    return records

def get_all_activities_for_compound(compound_id, 
    filter_pa_pi=True):
    print ("getting all activities for compound", compound_id)

    categories = get_categories()
    category_targets = {category: get_targets_for_category(category)
        for category in categories}

    conn = connect_to_mysqldb()


    category_activities = {category: 
        mysql_query(" UNION ALL ".join(
        (f"SELECT Pa, Pi from `{target_id}` where coconut_id='{compound_id}'" 
            for target_id in category_targets[category].values())), existing_conn=conn)
        for category in categories
    }

    conn.close()

    # category_activities = [("ALL",
    #     [(target, 
    #         activities[target_id][0], activities[target_id][1], #Pa, Pi
    #             activities[target_id][0] - activities[target_id][1]) # Pa-Pi
    #     for target, target_id in targets.items()
    #     if not filter_pa_pi or activities[target_id][0] > activities[target_id][1]]
    # )]

    print ("completed query")
    print ("building records")

    category_activities = [(category,
        [(target, target_activities[0], target_activities[1], #Pa, Pi
            target_activities[0] - target_activities[1]) # Pa-Pi
                for target, target_activities in zip(category_targets[category].keys(), category_activities_)
                    if not filter_pa_pi or target_activities[0] > target_activities[1]])
        for category, category_activities_ in category_activities.items()
    ]
    
    return category_activities
    

def get_compound_info(
    compound_id, 
    projection=None,
    get_activities=True,
    compound_info_collection="uniqueNaturalProduct",
    # pass_collection=PASS_COLLECTION
    ):

    print ("querying COCONUT database for info",
        "about compound with compound id", compound_id)

    db = connect_to_mongodb()

    coconut_collection = db[compound_info_collection]

    compound_info = coconut_collection.find_one(
        {"coconut_id": compound_id},
        projection=projection)

    if not get_activities:
        return compound_info

    # pass_collection = db[pass_collection]

    # categories = get_categories()

    # # perform query and keep all categories
    # compound_activities = pass_collection.find_one(
    #     {"coconut_id": compound_id}, 
    #     {"PASS_" + category: 1 for category in categories}.update({"PASS_ERROR": 1 })
    # )

    # if pd.isna(compound_activities["PASS_ERROR"]):
    #     pass_activities = []
    #     for category in categories:
    #         category_activities = compound_activities["PASS_" + category]

    #         pass_activities.append(
    #             (
    #                 category, 
    #                 [(target, activity["Pa"], activity["Pi"], activity["Pa"] - activity["Pi"])
    #                     for target, activity in category_activities.items()
    #                     if activity["Pa"] > activity["Pi"]]
    #             )
    #         )
    # else:
    #     pass_activities = None

    # db.client.close()
    pass_activities = get_all_activities_for_compound(compound_id)

    return compound_info, pass_activities

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

def write_smiles_to_file(
    username,
    records,
    static_dir="natural_products/static",
    output_dir="smiles",
    ):  

    output_dir = os.path.join(static_dir, output_dir, username)

    os.makedirs(output_dir, exist_ok=True)

    smiles_filename = os.path.join(output_dir,
        "results.smi")
    print ("writing smiles to", smiles_filename)

    # columns are ["coconut_id",  "molecule_name", "molecular_formula", "smiles", "Pa", "Pi", "Pa-Pi"]
    smiles = [(record[0], record[3]) for record in records]

    write_smiles(smiles, smiles_filename)

    return smiles_filename

def write_records_to_file(
    username,
    records,
    static_dir="natural_products/static",
    output_dir="results",
    ):

    output_dir = os.path.join(static_dir, output_dir, username)

    os.makedirs(output_dir, exist_ok=True)
    records = pd.DataFrame.from_records(records, 
        columns=["coconut_id",  "molecule_name", "molecular_formula", "smiles", "Pa", "Pi", "Pa-Pi"])

    records_filename = os.path.join(output_dir,
        "results.csv")
    print ("writing records to", records_filename)
    records.to_csv(records_filename)

    return records_filename
    

# def migrate_to_mysql():

#     db = connect_to_mongodb()
#     collection = db[PASS_COLLECTION]

#     target_map_inv = load_json("target_ids.json")

#     # create tables in MYSQL
#     for category in [
#         "EFFECTS",
#         "MECHANISMS",
#         'TOXICITY', 
#         'ANTITARGETS',
#         'METABOLISM', 
#         'GENE_EXPRESSION', 
#         'TRANSPORTERS'
#     ]:

#         print ("processing category", category)

#         targets = get_targets_for_category(category)

#         for i, target in enumerate(targets):
#             print ("processing target", target)

#             create_table = f"CREATE TABLE `{target_map_inv[target]}` (coconut_id VARCHAR(255) PRIMARY KEY, Pa SMALLINT, Pi SMALLINT )"
#             mysql_create_table(create_table)

#             # get existing records
#             existing_ids_query = f"select coconut_id from `{target_map_inv[target]}`"
#             existing_ids = mysql_query(existing_ids_query)
#             existing_ids = sorted({record[0] for record in existing_ids})

#             records = collection.find(
#                 {"$and": [{"coconut_id": {"$nin": existing_ids}},
#                     {f"PASS_{category}": {"$ne": np.NaN}}]},
#                 projection={"_id": 0, "coconut_id": 1, f"PASS_{category}.{target}": 1})

#             # insert records for that target
#             print (f"inserting records for target {target}")
#             sql = f"INSERT INTO `{str(target_map_inv[target])}`" + " (coconut_id, Pa, Pi) VALUES (%s, %s, %s)"
#             print ("using SQL command", sql)

#             print ("building rows to insert into SQL")
#             rows = [
#                 (record["coconut_id"], 
#                     int(record[f"PASS_{category}"][target]["Pa"]), 
#                     int(record[f"PASS_{category}"][target]["Pi"]), )
#                 for record in records
#             ]
#             mysql_insert_many(sql, rows)

#             print ("completed target", target, i+1, "/", len(targets))
#             print ("################")
#             print ()

#         print ("completed category", category)
#         print ("################")
#         print ()

if __name__ == "__main__":

    # hits = query_target_hits("Diarrhea", 500)

    # print (hits)

    compound_info, activities = get_compound_info("CNP0000002")
