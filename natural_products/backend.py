import os

import numpy as np
import pandas as pd

from datetime import datetime

import pymongo
from pymongo import MongoClient
from pymongo.errors import ServerSelectionTimeoutError

from rdkit import Chem
from rdkit.Chem.Draw import MolToImage

import urllib.parse as urlparse

import json

import mysql.connector
from mysql.connector.errors import ProgrammingError, IntegrityError


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

def connect_to_mongodb(host=HOST, port=PORT, db=DB, timeout=5000):
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

    db = connect_to_mongodb()
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

    db = connect_to_mongodb()
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

    db = connect_to_mongodb()

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

    db = connect_to_mongodb()

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

def _as_batch(cursor, batch_size=50):
    # iterate over something (pymongo cursor, generator, ...) by batch. 
    # Note: the last batch may contain less than batch_size elements.
    batch = []
    try:
        while True:
            for _ in range(batch_size):
                batch.append(next(cursor))
                if len(batch) % 1000 == 0:
                    print ("building batch...")
                    print ("current size:", len(batch))
            print ("yielding batch of size", len(batch))
            yield batch
            batch = []
    except StopIteration as e:
        if len(batch):
            yield batch

def mysql_query(query):

    mydb = mysql.connector.connect(
        host="192.168.0.49",
        user="david",
        password="c423612k",
        database="npaiengine"
    )

    mycursor = mydb.cursor()
    mycursor.execute(query)

    records = mycursor.fetchall()

    mydb.close()

    return records

def mysql_create_table(create_table):
    mydb = mysql.connector.connect(
        host="192.168.0.49",
        user="david",
        password="c423612k",
        database="npaiengine"
    )

    mycursor = mydb.cursor()
    try:
        mycursor.execute(create_table)
        print ("executed command", create_table)

    except ProgrammingError as e: # table already exists
        print ("table already exists")
        pass
    return 0

def mysql_insert_many(sql, rows):

    print ("inserting", len(rows), "rows")

    mydb = mysql.connector.connect(
        host="192.168.0.49",
        user="david",
        password="c423612k",
        database="npaiengine"
    )

    mycursor = mydb.cursor()
    mycursor.executemany(sql, rows)
    mydb.commit()

    mydb.close()

    return 0

def migrate_to_mysql():

    db = connect_to_mongodb()
    collection = db[PASS_COLLECTION]

    targets = sorted({target for category in get_categories() for target in get_targets_for_category(category)})
    target_map = {i: target for i, target in enumerate(targets)}
    target_map_inv = {target: i for i, target in target_map.items()}

    # create tables in MYSQL
    for category in get_categories()[::-1]:

        print ("processing category", category)

        targets = get_targets_for_category(category)

        for i, target in enumerate(targets):
            print ("processing target", target)

            create_table = f"CREATE TABLE `{target_map_inv[target]}` (coconut_id VARCHAR(255) PRIMARY KEY, Pa SMALLINT, Pi SMALLINT )"
            mysql_create_table(create_table)

            # get existing records
            existing_ids_query = f"select coconut_id from `{target_map_inv[target]}`"
            existing_ids = mysql_query(existing_ids_query)
            existing_ids = sorted({record[0] for record in existing_ids})

            records = collection.find(
                {"$and": [{"coconut_id": {"$nin": existing_ids}},
                    {f"PASS_{category}": {"$ne": np.NaN}}]},
                projection={"_id": 0, "coconut_id": 1, f"PASS_{category}.{target}": 1})

            # insert records for that target
            print (f"inserting records for target {target}")
            sql = f"INSERT INTO `{str(target_map_inv[target])}`" + " (coconut_id, Pa, Pi) VALUES (%s, %s, %s)"
            print ("using SQL command", sql)

            print ("building rows to insert into SQL")
            rows = [
                (record["coconut_id"], 
                    int(record[f"PASS_{category}"][target]["Pa"]), 
                    int(record[f"PASS_{category}"][target]["Pi"]), )
                for record in records
            ]
            mysql_insert_many(sql, rows)

            print ("completed target", target, i+1, "/", len(targets))
            print ("################")
            print ()

        print ("completed category", category)
        print ("################")
        print ()

if __name__ == "__main__":
    migrate_to_mysql()
    
    # threshold = 100

    # db = connect_to_mongodb()
    # collection = db[PASS_COLLECTION]

    # categories = get_categories()

    # category = categories[0]

    # print ("using category", category)

    # targets = list(get_targets_for_category(category))

    # with open("target_map.json", "r") as f:
    #     target_map = json.load(f)

    # target = targets[0]

    # target_id = target_map[target]

    # mongo_records = collection.find(
    #     {"$and":[
    #         {f"PASS_{category}.{target}.Pa": {"$gt": threshold}}, 
    #         {"$expr": {"$gt":[f"$PASS_{category}.{target}.Pa", f"$PASS_{category}.{target}.Pi"]}}
    #     ]},
    #     projection={"_id": 0, "coconut_id": 1})
    
    # from timeit import default_timer

    # start_time = default_timer()

    # mongo_records = sorted({record["coconut_id"] for record in mongo_records})

    # print ("mongodb", default_timer() - start_time)

    # mydb = mysql.connector.connect(
    #     host="192.168.0.49",
    #     user="david",
    #     password="c423612k",
    #     database="npaiengine"
    # )

    # mycursor = mydb.cursor()

    # mycursor.execute(f"select coconut_id from `{target_id}` where Pa>{threshold} and Pa>Pi")

    # start_time = default_timer()

    # sql_records = {record[0] for record in mycursor.fetchall()}

    # print ("SQL", default_timer() - start_time)

    # print (sorted(mongo_records)[:10])
    # print (sorted(sql_records)[:10])

    # print (len(mongo_records), len(sql_records))

    # assert len(mongo_records) == len(sql_records), "length"
    # # assert all([record in mongo_records for record in sql_records]), "sql in mongo"
    # # assert all([record in sql_records for record in mongo_records]), "mongo in sql"

    # print ("pass")

    # # record = collection.find_one({})#{"PASS_EFFECTS": {"$ne": np.NaN}})

    # # print (record)


    # category = "TOXICITY"

    # targets = get_targets_for_category(category)

    # print (len(targets))

    # for target in targets:
        # print (target)

    # from pymongo import IndexModel, ASCENDING, DESCENDING

    # # # target = "Yawning"
    # # # target = "PASS_{}.{}".format(category, target)
    # # # print ("processing target", target)
    # collection.create_indexes([
    #     IndexModel([(f"PASS_{category}.{target}.Pa", ASCENDING) for target in get_targets_for_category(category)],
    #         name=category)
    #         # for target in get_targets_for_category(category) #("Diarrhea", "Yawning", "Delusion")
    #     ], 
    # )
