import sys
import os.path
sys.path.insert(1, 
    os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))

import mysql.connector
from mysql.connector.errors import ProgrammingError, IntegrityError

from utils.io import load_json
from utils.pass_utils import remove_invalid_characters, parse_pass_spectra, get_all_targets, get_categories, get_targets_for_category

from functools import partial

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

def mysql_query(query, existing_conn=None):

    if existing_conn is None:
        db = connect_to_mysqldb()
    else:
        db = existing_conn

    mycursor = db.cursor()
    mycursor.execute(query)

    records = mycursor.fetchall()

    if existing_conn is None:
        db.close()

    return records

def mysql_create_table(create_table, existing_conn=None):
    if existing_conn is None:
        db = connect_to_mysqldb()
    else:
        db = existing_conn

    mycursor = db.cursor()
    try:
        mycursor.execute(create_table)
        print ("executed command", create_table)

    except ProgrammingError as e: # table already exists
        # print ("table already exists")
        print (e)
        pass

    if existing_conn is None:
        db.close()

    return 0

def mysql_insert_many(sql, rows, existing_conn=None):

    print ("inserting", len(rows), "rows")

    if existing_conn is None:
        db = connect_to_mysqldb()
    else:
        db = existing_conn

    mycursor = db.cursor()
    mycursor.executemany(sql, rows)
    mydb.commit()

    if existing_conn is None:
        db.close()

    return 0

def create_tables(): # too many columns

    targets = get_all_targets()
    n_targets = len(targets)

    for table in ("pa", "pi"):

        create_table = f"CREATE TABLE `{table}` (coconut_id VARCHAR(255) PRIMARY KEY"
        for target_id in range(n_targets):
            create_table += f", `{target_id}` SMALLINT"
        create_table += ")"

        mysql_create_table(create_table)

    return 0


def migrate_SDF_to_mysql(sdf_file):
    print ("reading SDF file from", sdf_file, "and inserting into SQL database")

    target_to_id = get_all_targets

    categories = get_categories()

    import pandas as pd

    # from rdkit.Chem.PandasTools import LoadSDF

    from collections import defaultdict

    from utils.rdkit_utils import LoadSDF

    chunks = LoadSDF(sdf_file, smilesName="SMILES", molColName=None, chunksize=5000)

    for chunk_no, chunk in enumerate(chunks):
       
        chunk = chunk.set_index("coconut_id", drop=True)
        chunk = chunk.loc[pd.isnull(chunk["PASS_ERROR"])] # find only valid compounds
        chunk = chunk[[f"PASS_{category}" for category in categories]] # drop unnecessary columns

        # create tables in MySQL
        for category in categories:

            print ("processing category", category)

            col = f"PASS_{category}"
            assert col in chunk.columns
            print ("parsing PASS activities")
            category_col = chunk[col].map(lambda s: 
                list(map(partial(parse_pass_spectra, mapping=target_to_id), 
                    map(remove_invalid_characters, s.split("\n")))),
                na_action="ignore")
            del chunk[col]

            # targets = get_targets_for_category(category)
            print ("building list of records for each target")
            targets = defaultdict(list)
            compound_ids = []
            for compound, target_activities in category_col.items():
                compound_ids.append(compound)
                for target, activities in target_activities:
                    targets[target].append((activities["Pa"], activities["Pi"]))
            del category_col

            conn = connect_to_mysqldb()

            for i, target_id in enumerate(targets):
                # assert target in target_to_id
                # target_id = target_to_id[target]
                print ("processing target_id", target_id)

                create_table = f"CREATE TABLE `{target_id}` (coconut_id VARCHAR(255) PRIMARY KEY, Pa SMALLINT, Pi SMALLINT )"
                mysql_create_table(create_table)

                # get existing records
                existing_ids_query = f"select coconut_id from `{target_id}`"
                existing_ids = mysql_query(existing_ids_query)
                existing_ids = {record[0] for record in existing_ids}

                # insert records for that target
                print (f"inserting records for target_id {target_id}")
                sql = f"INSERT INTO `{target_id}`" + " (coconut_id, Pa, Pi) VALUES (%s, %s, %s)"
                print ("using SQL command", sql)

                print ("building rows to insert into SQL")
                rows = [
                    (compound, *row) for compound, row in zip(compound, targets[target_id])
                    if compound not in existing_ids # compound_id
                ]
                mysql_insert_many(sql, rows, existing_conn=conn)

                print ("completed target_id", target_id, i+1, "/", len(targets))
                print ("################")
                print ()

            conn.close()

            print ("completed category", category)
            print ("################")
            print ()
       
        print ("completed chunk number", chunk_no)
        print ("################")
        print ()


if __name__ == "__main__":
    # mydb = connect_to_mysqldb()

    # print (mydb)

    # migrate_SDF_to_mysql("/media/david/26FE51D21F197BDF/coconut/COCONUT_0 (PASS2019).SDF")

    # create_tables()
    for category in get_categories():

        targets = get_targets_for_category(category)

        print (category, len(targets))