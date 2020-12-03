import sys
import os.path
sys.path.insert(1, 
    os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))

import pymongo
from pymongo import MongoClient
from pymongo.errors import ServerSelectionTimeoutError

from utils.io import load_json

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

def count_documents(collection, filter={}):
    return collection.count_documents(filter=filter)

def clear_collection(collection, filter={}):
    collection.remove(filter=filter)

if __name__ == "__main__":
    db = connect_to_mongodb()
    print (db)