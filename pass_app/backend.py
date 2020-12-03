
import os
import sys
import os.path
sys.path.insert(1, 
    os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))

from django.utils.encoding import smart_str

from datetime import datetime

from pymongo import MongoClient

import numpy as np
import pandas as pd

from rdkit.Chem.PandasTools import LoadSDF

from utils.email_utils import send_mail
from utils.mongodb_utils import connect_to_mongodb

# HOST = "192.168.0.49"
# PORT = 27017
# DB = "PASS_TEST"

# def connect_to_mongodb()():
    # print ("connecting to MongoDB database", 
        # DB, "using host",
        # HOST, "and port", PORT)
    # client = MongoClient(HOST, PORT)
    # return client[DB]

def remove_invalid_characters(s):
    return s.replace(".", "")

def parse_pass_spectra(s, ):
    split = s.split()
    pa = split[0]
    pi = split[1]
    activity = " ".join(split[2:])
    return {activity: {"Pa": pa, "Pi": pi}}

def write_PASS_hits_to_db(
    pass_file,
    collection="PASS_hits"):
    assert isinstance(pass_file, str)
    assert pass_file.endswith(".sdf")
    print ("reading PASS hits from", pass_file,
        "and extracting activity data")

    pass_activities = LoadSDF(pass_file, 
        smilesName="SMILES", molColName=None)

    db = connect_to_mongodb()()

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
                value = list(
                    map(parse_pass_spectra, 
                        map(remove_invalid_characters,
                            row[col].split("\n"))))
            else:
                value = row[col]
            entry.update({col: value})
        entry.update({"time": str(datetime.now())})
        records.append(entry)

    print ("inserting", len(records), "records")
    hit_collection.insert_many(records)

    db.client.close()

def pass_predict(
    user_name,
    user_email,
    sdf_file, 
    output_dir=os.path.join("pass_app", "static", "pass_app",
        "files")):

    output_dir = os.path.join(output_dir, user_email)
    os.makedirs(output_dir, exist_ok=True)

    filename = sdf_file.name
    assert filename.endswith(".sdf")
    filename = os.path.splitext(filename)[0]

    output_dir = os.path.join(output_dir, filename)
    os.makedirs(output_dir, exist_ok=True)

    input_file = os.path.join(output_dir,
        filename + "-in.sdf")
    with open(input_file, "wb+") as out_file:
        for chunk in sdf_file.chunks():
            out_file.write(chunk)

    pass_out_file = os.path.join(output_dir, 
        filename +"-PASS-out.sdf")

    cmd = "PASS2019toSDF.exe {} {}".format(input_file, pass_out_file)
    print ("executing command:", cmd)

    ret = os.system(cmd)
    assert ret == 0
    assert os.path.exists(pass_out_file)

    # write PASS spectra to database 
    write_PASS_hits_to_db(pass_out_file)

    print ("removing", input_file)
    os.remove(input_file)

    # return smart_str(pass_out_file)

    # send mail containing results

    send_mail(user_name, user_email, pass_out_file)

    return 0

if __name__ == "__main__":
    name = "david"
    email = "davemcdonald93@gmail.com"
    sdf_file = "/home/david/Desktop/test.sdf"



