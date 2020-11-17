
import os

from django.utils.encoding import smart_str

from datetime import datetime

from pymongo import MongoClient

import numpy as np
import pandas as pd

from rdkit.Chem.PandasTools import LoadSDF

HOST = "192.168.0.49"
PORT = 27017
DB = "PASS_TEST"

def connect_to_db():
    print ("connecting to MongoDB database", 
        DB, "using host",
        HOST, "and port", PORT)
    client = MongoClient(HOST, PORT)
    return client[DB]

def count_documents(collection, filter={}):
    return collection.count_documents(filter=filter)

def clear_collection(collection, filter={}):
    collection.remove(filter=filter)

def parse_pass_spectra(s, ):
    split = s.split()
    return {"Pa": split[0], "Pi": split[1], "activity" : " ".join(split[2:])}

def write_PASS_hits_to_db(
    pass_file,
    collection="PASS_hits"):
    assert isinstance(pass_file, str)
    assert pass_file.endswith(".sdf")

    pass_activities = LoadSDF(pass_file, 
            smilesName='SMILES', molColName=None)

    db = connect_to_db()

    print ("writing PASS hits to MongoDB:", DB, 
        "collection:", collection)

    hit_collection = db[collection]

    records = []

    print ("determining hits")

    for _, row in pass_activities.iterrows():
        
        entry = {}
        for col in row.index:
            if pd.isnull(row[col]):
                continue
            if col == "PASS_ACTIVITY_SPECTRUM":
                value = list(map(parse_pass_spectra, row[col].split("\n")))
            else:
                value = row[col]
            entry.update({col: value})
        entry.update({"time": str(datetime.now()))
        records.append(entry)

    print ("inserting", len(records), "records")
    hit_collection.insert_many(records)

    db.client.close()

def handle_uploaded_file(f):

    filename = f.name
    assert filename.endswith(".sdf")

    temp_file = "temp.sdf"
    with open(temp_file, "wb+") as out_file:
        for chunk in f.chunks():
            out_file.write(chunk)

    pass_out_file = os.path.splitext(filename)[0] +\
        "-PASS-out.sdf"

    cmd = "PASS2019toSDF {} {}".format(temp_file, pass_out_file)
    print ("executing command:", cmd)

    ret = os.system(cmd)

    assert ret == 0

    # write PASS spectra to database 
    write_PASS_hits_to_db(pass_out_file)

    return smart_str(pass_out_file)



