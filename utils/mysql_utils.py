import sys
import os.path
sys.path.insert(1, 
    os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))

import mysql.connector
from mysql.connector.errors import ProgrammingError, IntegrityError

from utils.io import load_json
from utils.pass_utils import remove_invalid_characters, parse_pass_spectra, get_all_targets, get_categories, get_targets_for_category, get_all_compounds

from functools import partial

import re

def load_mysql_credentials():
    return load_json("mysql_credentials.json")

def connect_to_mysqldb(host=None, user=None, password=None, database=None):
    mysql_credentials = load_mysql_credentials()
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

def mysql_execute(cmd, existing_conn=None):


    if existing_conn is None:
        db = connect_to_mysqldb()
    else:
        db = existing_conn

    print ("executing MySQL cmd:", cmd)

    cursor = db.cursor()
    cursor.execute(cmd)

    if existing_conn is None:
        db.close()

    return 0

def mysql_query(query, return_cols=False, existing_conn=None):

    if existing_conn is None:
        db = connect_to_mysqldb()
    else:
        db = existing_conn

    print ("executing MySQL query:", query)

    cursor = db.cursor()
    cursor.execute(query)

    records = cursor.fetchall()
    print ("NUMBER OF HITS:", len(records))

    if return_cols:
        records = records, cursor.column_names

    if existing_conn is None:
        db.close()

    print ()

    return records

def mysql_create_table(create_table, existing_conn=None):
    if existing_conn is None:
        db = connect_to_mysqldb()
    else:
        db = existing_conn

    cursor = db.cursor()
    try:
        cursor.execute(create_table)
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

    cursor = db.cursor()

    if isinstance(rows, list):
        print ("inserting", len(rows), "row(s)")
        cursor.executemany(sql, rows)
        db.commit()

    else:

        for chunk in to_chunks(rows):
            print ("inserting", len(chunk), "row(s)")
            cursor.executemany(sql, chunk)
            db.commit()


    if existing_conn is None:
        db.close()

    return 0

def sanitise_names(names):
    return  [
        re.sub(r"( |\+|-|\*|/|=|<|>|\(|\)|,|\.|'|\[|\]|:|;|{|})", "_", name)
        for name in names]

if __name__ == "__main__":
    
    query = '''
    SELECT acc
    FROM uniprot
    LIMIT 10000

    ''' 

    records = mysql_query(query)

    with open("uniprots.txt", "w") as f:
        f.write("\n".join((record[0] for record in records)))