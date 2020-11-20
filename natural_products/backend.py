import os

from datetime import datetime

import pymongo
from pymongo import MongoClient

from rdkit import Chem
from rdkit.Chem.Draw import MolToImage

HOST = "192.168.0.49"
PORT = 27017
DB = "COCONUT"

def connect_to_db(host=HOST, port=PORT, db=DB):
    print ("connecting to MongoDB database", 
        db, "using host",
        host, "and port", port)
    client = MongoClient(host, port)
    return client[db]

def query_pass(
    target, 
    threshold,
    collection="PASS"):

    print ("querying PASS database for compounds",
        "that hit target", target, "greater than or equal to",
        "threshold", threshold)

    db = connect_to_db()
    pass_collection = db[collection]

    query = {target: {"$gte": threshold}, }
    filter_ = {"_id": 0, "coconut_id": 1, 
        "molecular_formula": 1, target: 1}

    cursor = pass_collection.find(query, filter_)\
        .sort([(target, pymongo.DESCENDING)]) # descending
    print ("performed query")

    print ("iterating over query")

    records = [
        (record["coconut_id"], 
            record["molecular_formula"], record[target])
        for record in cursor
    ]

    db.client.close()

    return records

def get_compound_info(
    compound_id, ):

    print ("querying COCONUT database for info",
        "about compound with compound id", compound_id)

    db = connect_to_db()

    coconut_collection = db["uniqueNaturalProduct"]

    compound_info = coconut_collection.find_one(
        {"coconut_id": compound_id})

    pass_collection = db["PASS"]
    pass_activities = pass_collection.find_one(
        {"coconut_id": compound_id})
    
    db.client.close()

    assert compound_info is not None, compound_id
    assert pass_activities is not None

    return compound_info, pass_activities

def draw_molecule(smiles, 
        img_filename="natural_products/temp.png", 
        static_dir="natural_products/static"):
    mol = Chem.MolFromSmiles(smiles)
    img = MolToImage(mol)

    img.save(os.path.join(static_dir,
        img_filename))

    return img_filename
