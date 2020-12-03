import sys
import os.path
sys.path.insert(1, 
    os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))

import mysql.connector
from mysql.connector.errors import ProgrammingError, IntegrityError

from utils.io import load_json
from utils.mongodb_utils import connect_to_mongodb

def connect_to_mysqldb(host=None, user=None, password=None, database=None):
    mysql_credentials = load_json("mysql_credentials.json")
    if host is None:
        host = mysql_credentials["host"]
    if user is None:
        user = mysql_credentials["user"]
    if password is None:
        password = mysql_credentials["password"]
    if database is None:
        database = mysql_credentials["database"]
    print ("connecting to MySQL database", 
        database, "using host",
        host, "and user", user)
    mydb = mysql.connector.connect(
        host=host,
        user=user,
        password=password,
        database=database
    )
    return mydb

def mysql_query(query):

    mydb = connect_to_mysqldb()

    mycursor = mydb.cursor()
    mycursor.execute(query)

    records = mycursor.fetchall()

    mydb.close()

    return records

def mysql_create_table(create_table):
    mydb = connect_to_mysqldb()

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

    mydb = connect_to_mysqldb()

    mycursor = mydb.cursor()
    mycursor.executemany(sql, rows)
    mydb.commit()

    mydb.close()

    return 0

# def migrate_to_mysql(PASS_COLLECTION="PASS_ALL"):

#     db = connect_to_mongodb()
#     collection = db[PASS_COLLECTION]

#     target_map_inv = load_json("target_ids.json")

#     # create tables in MySQL
#     for category in [
#         # "EFFECTS",
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
    mydb = connect_to_mysqldb()

    # print (mydb)

    # migrate_to_mysql()