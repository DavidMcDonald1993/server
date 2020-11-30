import os

import pandas as pd

from datetime import datetime

import pymongo
from pymongo import MongoClient
from pymongo.errors import ServerSelectionTimeoutError

from rdkit import Chem
from rdkit.Chem.Draw import MolToImage

import urllib.parse as urlparse

HOST = "192.168.0.49"
PORT = 27017
DB = "COCONUT"
PASS_COLLECTION = "PASS_ALL"

target_categories = [
    "EFFECTS",
    "MECHANISMS",
    'TOXICITY', 'ANTITARGETS',
    'METABOLISM', 'GENE_EXPRESSION', 'TRANSPORTERS'
]

def connect_to_db(host=HOST, port=PORT, db=DB, timeout=5000):
    print ("connecting to MongoDB database", 
        db, "using host",
        host, "and port", port)
    client = MongoClient(host, port, serverSelectionTimeoutMS=timeout)
    return client[db]

# def query_database(collection_name, query, show, return_one):

#     for host in ("localhost", HOST):
#         try:
#             db = connect_to_db(host=host)
#             collection = db[collection_name]

#             if return_one:
#                 query_fun = collection.find_one
#             else:
#                 query_fun = collection.find 

#             print ("performing query with filter:", query,
#                 "and showing", show)

#             records = query_fun(query, show)

#             db.client.close()

#             return records

#         except ServerSelectionTimeoutError as e:
#             pass

#     raise ServerSelectionTimeoutError

def get_categories():
    return target_categories

def get_targets_for_category(category,
    collection=PASS_COLLECTION):
    assert category in target_categories

    print ("identifying targets for category", category)

    db = connect_to_db()
    pass_collection = db[collection]

    query = {"coconut_id": "CNP0000002"}
    show = {"_id": 0, "PASS_" + category: 1}

    record = pass_collection.find_one(query, show)

    # record = query_database(collection_name=collection,
        # query=query, show=show, return_one=True)
  
    return record["PASS_" + category].keys()


def query_pass_activities(
    category,
    target, 
    threshold,
    filter_pa_pi=True,
    collection=PASS_COLLECTION):
    assert category in target_categories
    category = "PASS_" + category

    print ("querying PASS database for compounds",
        "that hit target", target, "in category", category,
         "with Pa greater than or equal to",
        "threshold", threshold)

    db = connect_to_db()
    pass_collection = db[collection]

    category_target = "{}.{}".format(category, target)

    query = [{"{}.Pa".format(category_target): {"$gte": threshold}, }]
    if filter_pa_pi:
        query.append( 
            {"$expr": { "$gt": 
                ["${}.Pa".format(category_target), "${}.Pi".format(category_target)] } } )
    
    show = {
        "_id": 0, 
        "coconut_id": 1, 
        "name": 1,
        "molecular_formula": 1, 
        # "SMILES": 1,
        category_target: 1
    }

    query = {"$and": query}

    print ("filtering with", query)
    print ("showing", show)

    cursor = pass_collection.find(query, show)#\
        # .sort([(target+".Pa", pymongo.DESCENDING), ]) 

    # cursor = query_database(collection_name=collection,
        # query=query, show=show, return_one=False)
    
    cursor.batch_size(1000000)

    print ("iterating over records")

    records = [
        (record["coconut_id"], 
            (urlparse.unquote(record["name"]).capitalize() 
                if record["name"] is not None else "<NO NAME>"),
            record["molecular_formula"], 
            # record["SMILES"], 
            record[category][target]["Pa"], record[category][target]["Pi"], 
            record[category][target]["Pa"] - record[category][target]["Pi"])
        for record in cursor
    ]

    db.client.close()

    print ("built records", len(records))

    return records


def get_all_compounds(    
    compound_info_collection="uniqueNaturalProduct",):

    print ("querying COCONUT database for info",
        "about all compounds")

    db = connect_to_db()

    coconut_collection = db[compound_info_collection]

    query = {}
    show = {
        "_id": 0, 
        "coconut_id": 1, 
        "name": 1,
        "molecular_formula": 1, 
        "clean_smiles": 1
    }

    cursor = coconut_collection.find(query, show)\
        .sort([("coconut_id", pymongo.ASCENDING)])
    cursor.batch_size(1000000)

    print ("iterating over query")

    records = [
        (record["coconut_id"], 
            (urlparse.unquote(record["name"]).capitalize() 
                    if "name" in record else "<NO NAME>"),
            record["molecular_formula"],
            record["clean_smiles"])
        for record in cursor
    ]

    db.client.close()

    print ("built records", len(records))

    return records



def get_compound_info(
    compound_id, 
    compound_info_collection="uniqueNaturalProduct",
    pass_collection=PASS_COLLECTION):

    print ("querying COCONUT database for info",
        "about compound with compound id", compound_id)

    db = connect_to_db()

    coconut_collection = db[compound_info_collection]

    compound_info = coconut_collection.find_one(
        {"coconut_id": compound_id})

    pass_collection = db[pass_collection]

    categories = get_categories()

    # perform query and keep all categories
    compound_activities = pass_collection.find_one(
        {"coconut_id": compound_id}, 
        {"PASS_" + category: 1 for category in categories}.update({"PASS_ERROR": 1 })
    )

    if pd.isna(compound_activities["PASS_ERROR"]):
        pass_activities = []
        for category in categories:
            category_activities = compound_activities["PASS_" + category]

            pass_activities.append(
                (
                    category, 
                    [(target, activity["Pa"], activity["Pi"], activity["Pa"] - activity["Pi"])
                        for target, activity in category_activities.items()
                        if activity["Pa"] > activity["Pi"]]
                )
            )
    else:
        pass_activities = None

    db.client.close()

    assert compound_info is not None, compound_id

    return compound_info, pass_activities

if __name__ == "__main__":

    compound_info, activities = get_compound_info("CNP0000002")

    print (activities[0])

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
