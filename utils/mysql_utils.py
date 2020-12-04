import sys
import os.path
sys.path.insert(1, 
    os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))

import mysql.connector
from mysql.connector.errors import ProgrammingError, IntegrityError

from utils.io import load_json
# from utils.mongodb_utils import connect_to_mongodb
from utils.pass_utils import remove_invalid_characters, parse_pass_spectra

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

def migrate_SDF_to_mysql(sdf_file):
    print ("reading SDF file from", sdf_file, "and inserting into SQL database")

    target_to_id = load_json("target_ids.json")

    categories = [
        "EFFECTS",
        "MECHANISMS",
        'TOXICITY', 
        'ANTITARGETS',
        'METABOLISM', 
        'GENE_EXPRESSION', 
        'TRANSPORTERS',
    ]

    import pandas as pd

    # from rdkit.Chem.PandasTools import LoadSDF

    from collections import defaultdict

    from utils.rdkit_utils import LoadSDF

    chunks = LoadSDF(sdf_file, smilesName="SMILES", molColName=None, chunksize=5000)

    for chunk_no, chunk in chunks:
       
        chunk = chunk.set_index("coconut_id", drop=True)
        chunk = chunk.loc[pd.isnull(chunk["PASS_ERROR"])] # find only valid compounds
        chunk = chunk[[f"PASS_{category}" for category in categories]] # drop unnecessary columns

        # create tables in MySQL
        for category in categories:

            print ("processing category", category)

            col = f"PASS_{category}"
            assert col in chunk.columns
            category_col = chunk[col].map(lambda s: list(map(parse_pass_spectra, 
                            map(remove_invalid_characters, s.split("\n")))),
                na_action="ignore")

            # targets = get_targets_for_category(category)
            targets = defaultdict(list)
            for compound, target_activities in category_col.items():
                for target, activities in target_activities:
                    targets[target].append((compound, activities["Pa"], activities["Pi"]))

            for i, target in enumerate(targets):
                assert target in target_to_id
                target_id = target_to_id[target]
                print ("processing target", target)

                create_table = f"CREATE TABLE `{target_id}` (coconut_id VARCHAR(255) PRIMARY KEY, Pa SMALLINT, Pi SMALLINT )"
                mysql_create_table(create_table)

                # get existing records
                existing_ids_query = f"select coconut_id from `{target_id}`"
                existing_ids = mysql_query(existing_ids_query)
                existing_ids = {record[0] for record in existing_ids}

                # insert records for that target
                print (f"inserting records for target {target}")
                sql = f"INSERT INTO `{target_id}`" + " (coconut_id, Pa, Pi) VALUES (%s, %s, %s)"
                print ("using SQL command", sql)

                print ("building rows to insert into SQL")
                rows = [
                    row for row in targets[target]
                    if row[0] not in existing_ids # compound_id
                ]
                mysql_insert_many(sql, rows)

                print ("completed target", target, i+1, "/", len(targets))
                print ("################")
                print ()

            print ("completed category", category)
            print ("################")
            print ()
       
        print ("completed chunk number", chunk_no)
        print ("################")
        print ()


if __name__ == "__main__":
    # mydb = connect_to_mysqldb()

    # print (mydb)

    migrate_SDF_to_mysql("/media/david/26FE51D21F197BDF/coconut/COCONUT_0 (PASS2019).SDF")