import sys
import os.path
sys.path.insert(1, 
    os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))

import pymongo
from pymongo import MongoClient
from pymongo.errors import ServerSelectionTimeoutError

from utils.io import load_json, write_json

def connect_to_mongodb(host=None, port=None, db=None, timeout=5000):
    mongodb_credentials = load_json("mongodb_credentials.json")
    if host is None:
        host = mongodb_credentials["host"]
    if port is None:
        port = mongodb_credentials["port"]
    if db is None:
        db = mongodb_credentials["db"]
    print ("connecting to MongoDB database", 
        db, "using host",
        host, "and port", port)
    client = MongoClient(host, port, serverSelectionTimeoutMS=timeout)
    return client[db]

def count_documents(collection, filter_={}):
    return collection.count_documents(filter=filter_)

def clear_collection(collection, filter={}):
    collection.remove(filter=filter)

if __name__ == "__main__":

    sql = f'''
    SELECT compound_id, coconut_id
    FROM compounds 
    '''
    from utils.mysql_utils import mysql_query
    compound_to_id = {c: i for i, c in mysql_query(sql)}

    db = connect_to_mongodb()
    collection = db["uniqueNaturalProduct"]

    from collections import defaultdict
    from itertools import count

    import re

    cursor = collection.find({}, {"_id": 0, "coconut_id": 1, "textTaxa":1 })
    # n_records = len(cursor)
    # print ("retrieved", n_records, "from MongoDB")

    kingdom_to_id = defaultdict(count(1).__next__)
    species_to_id = defaultdict(count(1).__next__)

    compound_to_kingdom = defaultdict(list)
    compound_to_species = defaultdict(list)
    # kingdom_to_species = defaultdict(list)

    import unicodedata

    valid_kingdoms = {"marine", "plants", "fungi", "bacteria", "animals"}

    for i, record in enumerate(cursor):
        if i % 1000 == 0: print (i)
        coconut_id = record["coconut_id"]
        compound_id = compound_to_id[coconut_id]
        tax = record["textTaxa"]

        # print (tax)

        kingdoms = []
        species = []

        for t in tax:
            t = unicodedata.normalize('NFKD', t).encode('ascii','ignore').decode("utf-8")
            if t == "notax": continue
            # if " " not in t and re.match(r"^[a-z]", t):
            if t in valid_kingdoms:
                t = t.title()
                kingdoms.append(t)
            else:
                t = t.title()
                species.append(t)

        for k in kingdoms:
            kingdom_id = kingdom_to_id[k]
            compound_to_kingdom[compound_id].append(kingdom_id)
        for s in species:
            species_id = species_to_id[s]
            compound_to_species[compound_id].append(species_id)

    write_json(kingdom_to_id, "kingdom_to_id.json")
    write_json(species_to_id, "species_to_id.json")

    write_json(compound_to_kingdom, "compound_to_kingdom.json")
    write_json(compound_to_species, "compound_to_species.json")

    # unique_nps = sorted((record["coconut_id"] 
    #     for record in cursor))

    # unique_nps = {c: i for i, c in enumerate(unique_nps)}

    # write_json(unique_nps, "compound_ids.json")