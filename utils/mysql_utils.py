import sys
import os.path
sys.path.insert(1, 
    os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))

import mysql.connector
from mysql.connector.errors import ProgrammingError, IntegrityError

from utils.io import load_json
from utils.pass_utils import remove_invalid_characters, parse_pass_spectra, get_all_targets, get_categories, get_targets_for_category, get_all_compounds

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
    db = mysql.connector.connect(
        host=host,
        user=user,
        password=password,
        database=database,
        autocommit=False,
    )
    return db

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
        print (e)
        pass

    if existing_conn is None:
        db.close()

    return 0

def mysql_insert_many(sql, rows, existing_conn=None, chunksize=1000000):

    def to_chunks(rows):
        chunk = []

        for row in rows:
            chunk.append(row)
            if len(chunk) == chunksize:
                yield chunk
                chunk = []
        if len(chunk) > 0:
            yield chunk


    if existing_conn is None:
        db = connect_to_mysqldb()
    else:
        db = existing_conn

    mycursor = db.cursor()

    if isinstance(rows, list):
        print ("inserting", len(rows), "rows")
        mycursor.executemany(sql, rows)
    else:

        for chunk in to_chunks(rows):
            print ("inserting", len(chunk), "rows")
            mycursor.executemany(sql, chunk)

    db.commit()

    if existing_conn is None:
        db.close()

    return 0

def create_tables(): 

    compounds = get_all_compounds()
    targets = get_all_targets()
    # n_targets = len(targets)

    conn = connect_to_mysqldb()

    # create compounds table
    create_compounds_table = f"CREATE TABLE `compounds` (compound_id INT, coconut_id VARCHAR(255), PRIMARY KEY (compound_id))"
    mysql_create_table(create_compounds_table, existing_conn=conn)

    # insert into compounds table
    insert_compounds_sql = "INSERT INTO `compounds`" + " (compound_id, coconut_id) VALUES (%s, %s) ON DUPLICATE KEY UPDATE compound_id=compound_id"
    rows = [(i, c) for c, i in compounds.items()]
    mysql_insert_many(insert_compounds_sql, rows, existing_conn=conn)

    # create targets table
    create_targets_table = f"CREATE TABLE `targets` (target_id INT, target_name VARCHAR(255), PRIMARY KEY (target_id))"
    mysql_create_table(create_targets_table, existing_conn=conn)

    # insert into targets table
    insert_targets_sql = "INSERT INTO `targets`" + " (target_id, target_name) VALUES (%s, %s) ON DUPLICATE KEY UPDATE target_id=target_id"
    rows = [(i, t) for t, i in targets.items()]
    mysql_insert_many(insert_targets_sql, rows, existing_conn=conn)

    # create activities table
    create_activities_table = f"CREATE TABLE `activities` (compound_id INT, target_id INT, Pa SMALLINT, Pi SMALLINT, " +\
        "PRIMARY KEY (compound_id, target_id), FOREIGN KEY (compound_id) REFERENCES compounds(compound_id), FOREIGN KEY (target_id) REFERENCES targets(target_id))"
    mysql_create_table(create_activities_table, existing_conn=conn)

    # for table in ("pa", "pi"):

    #     create_table = f"CREATE TABLE `{table}` (coconut_id VARCHAR(255) PRIMARY KEY"
    #     for target_id in range(n_targets):
    #         create_table += f", `{target_id}` SMALLINT"
    #     create_table += ")"

    #     mysql_create_table(create_table)

    return 0


def migrate_SDF_to_mysql(sdf_file):
    print ("reading SDF file from", sdf_file, "and inserting into SQL database")

    compound_to_id = get_all_compounds()
    target_to_id = get_all_targets()

    categories = get_categories()

    import pandas as pd

    # from rdkit.Chem.PandasTools import LoadSDF

    from collections import defaultdict

    from utils.rdkit_utils import LoadSDF

    chunks = LoadSDF(sdf_file, smilesName="SMILES", molColName=None, chunksize=5000)

    sql = "INSERT INTO `activities`" + " (compound_id, target_id, Pa, Pi) VALUES (%s, %s, %s, %s) "+\
        "ON DUPLICATE KEY UPDATE compound_id=compound_id, target_id=target_id"

    for chunk_no, chunk in enumerate(chunks):

        chunk = chunk.loc[chunk["coconut_id"].isin(compound_to_id)]
       
        chunk = chunk.set_index("coconut_id", drop=True)
        if "PASS_ERROR" in chunk.columns:
            chunk = chunk.loc[pd.isnull(chunk["PASS_ERROR"])] # find only valid compounds
        chunk = chunk[[f"PASS_{category}" for category in categories]] # drop unnecessary columns
        assert not any(chunk.isnull().any(axis=0))

        print ("processing chunk", chunk_no)
        print ("chunk size:", chunk.shape[0])

        # create tables in MySQL
        for category in categories:

            print ("processing category", category)

            col = f"PASS_{category}"
            assert col in chunk.columns
            print ("parsing PASS activities")
            category_col = chunk[col].map(lambda s: 
                map(partial(parse_pass_spectra, mapping=target_to_id), 
                    map(remove_invalid_characters, s.split("\n"))),
                na_action="ignore")
            del chunk[col]

            # targets = get_targets_for_category(category)
            print ("building list of records for each target")
            targets = defaultdict(list)
            compound_ids = []
            for compound, target_activities in category_col.items():
                compound_ids.append(compound_to_id[compound])
                for target, activities in target_activities:
                    targets[target].append((activities["Pa"], activities["Pi"]))
            del category_col

            conn = connect_to_mysqldb()

            rows = (
                (compound, target_id, *row) 
                for target_id in targets
                for compound, row in zip(compound_ids, targets[target_id])
            )

            mysql_insert_many(sql, rows, existing_conn=conn)

            # for i, target_id in enumerate(targets):
            #     # assert target in target_to_id
            #     # target_id = target_to_id[target]
            #     # print ("processing target_id", target_id)

            #     # create_table = f"CREATE TABLE `{target_id}` (coconut_id VARCHAR(255) PRIMARY KEY, Pa SMALLINT, Pi SMALLINT )"
            #     # mysql_create_table(create_table, existing_conn=conn)

            #     # # get existing records
            #     # existing_ids_query = f"select coconut_id from `{target_id}`"
            #     # existing_ids = mysql_query(existing_ids_query, existing_conn=conn)
            #     # existing_ids = {record[0] for record in existing_ids}

            #     # insert records for that target
            #     print (f"inserting records for target_id {target_id}")
            #     sql = f"INSERT INTO `{target_id}`" + " (coconut_id, Pa, Pi) VALUES (%s, %s, %s) ON DUPLICATE KEY UPDATE coconut_id=coconut_id"
            #     print ("using SQL command", sql)

            #     print ("building rows to insert into SQL")
            #     rows = [
            #         (compound, *row) for compound, row in zip(compound_ids, targets[target_id])
            #         # if compound not in existing_ids # compound_id
            #     ]
            #     mysql_insert_many(sql, rows, existing_conn=conn)

            #     print ("completed target_id", target_id, i+1, "/", len(targets))
            #     print ("################")
            #     print ()

            conn.close()

            print ("completed category", category, "for chunk no", chunk_no)
            print ("################")
            print ()
       
        print ("completed chunk number", chunk_no)
        print ("################")
        print ()


if __name__ == "__main__":

    # create_tables()

    # conn = connect_to_mysqldb(database="test")

    # sql = f"INSERT INTO `users`" + " (firstname, surname) VALUES (%s, %s) ON DUPLICATE KEY UPDATE firstname = firstname"

    # rows = [
    #     ("Dora", "Emese"),
    #     ("Jack", "Statham"),
    #     ("Rebecca", "Bennett"),
    # ]

    # mysql_insert_many(sql, rows, existing_conn=conn)

    # records = mysql_query("SELECT * from users", existing_conn=conn)

    # print (records)

    migrate_SDF_to_mysql("/mnt/e/coconut/out/COCONUT_0 (PASS2019).SDF")

    # create_tables()

    # from timeit import default_timer

    # compound_id = "CNP0000002"

    # for category in ["EFFECTS"]:


    #     targets = get_targets_for_category(category)
    #     ids = list(targets.values())
    #     # anchor = ids[0]

    #     query = " UNION ALL ".join((f"SELECT Pa, Pi from `{i}` where coconut_id='{compound_id}'" for i in ids))
        # print (query)

        # # print (category, len(targets))

        # query = f"SELECT {anchor}.Pa, {anchor}.Pi"
        # for i in ids[1:]:
        #     query += f", {i}.Pa, {i}.Pi"
        # query += f" FROM `{anchor}`"
        # for i in ids[1:]:
        #     query += f", `{i}`"
        # query += f" WHERE `{i}`.coconut_id={compound_id} AND "
        # query += " AND ".join((f"`{anchor}`.coconut_id=`{i}`.coconut_id" for i in ids[1:]))


        # conn = connect_to_mysqldb()

        # start = default_timer()

        # records = mysql_query(query, conn)

        # print (default_timer() - start)

        # start = default_timer()

        # validation = [mysql_query(f"select Pa, Pi from `{i}` where coconut_id='{compound_id}'", conn)[0] for i in ids]

        # print (default_timer() - start)

        # # print (len(targets), len(ids), len(records), sum(validation))
        # assert len(records) == len(validation)

        # for r, v in zip(records, validation):
        #     # print (r, v)
        #     assert r[0]==v[0]
        #     assert r[1]==v[1]



        # break