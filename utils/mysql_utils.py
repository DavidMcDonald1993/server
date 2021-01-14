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

def mysql_query(query, return_cols=False, existing_conn=None):

    if existing_conn is None:
        db = connect_to_mysqldb()
    else:
        db = existing_conn

    print ("executing MySQL query:", query)

    cursor = db.cursor()
    cursor.execute(query)

    records = cursor.fetchall()

    if return_cols:
        records = records, cursor.column_names

    if existing_conn is None:
        db.close()

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
        print ("inserting", len(rows), "rows")
        cursor.executemany(sql, rows)
    else:

        for chunk in to_chunks(rows):
            print ("inserting", len(chunk), "rows")
            cursor.executemany(sql, chunk)

    db.commit()

    if existing_conn is None:
        db.close()

    return 0

def sanitise_names(names):
    return  [re.sub(r"( |\+|-|\*|/|=|<|>|\(|\)|,|\.|'|\[|\]|:|;)", "_", name)
        for name in names]

def get_all_targets_for_categories(categories=None, existing_conn=None):

    if categories is not None:
        if not isinstance(categories, str):
            if isinstance(categories, list) or isinstance(categories, set):
                categories = tuple(categories)
            assert isinstance(categories, tuple)
            if len(categories) == 1:
                categories = categories[0]

    query = f'''
       SELECT c.category_name, t.target_name
       FROM categories AS c
       INNER JOIN category_members AS m ON (c.category_id=m.category_id) 
       INNER JOIN targets AS t ON (m.target_id=t.target_id)
       {f"WHERE c.category_name IN {categories}" 
        if isinstance(categories, tuple) 
            else f"WHERE c.category_name='{categories}'" if isinstance(categories, str)
            else ""}
    '''

    return mysql_query(query, existing_conn=existing_conn)

def get_uniprots_for_targets(targets, existing_conn=None):
    if not isinstance(targets, str):
        if isinstance(targets, list) or isinstance(targets, set):
            targets = tuple(targets)
        assert isinstance(targets, tuple)
        print ("querying mySQL database for UNIPROTS associated with",
            len(targets), "targets")
        if len(targets) == 1:
            targets = targets[0]
    query = f'''
        SELECT DISTINCT t.target_name, u.acc, tu.score
        FROM targets AS t
        INNER JOIN targets_to_uniprot AS tu ON (t.target_id=tu.target_id)
        INNER JOIN uniprot AS u ON (tu.uniprot_id=u.uniprot_id)
        WHERE t.target_name {f"IN {targets}" if isinstance(targets, tuple)
            else f"='{targets}'"}
    '''
    return mysql_query(query, existing_conn=existing_conn)

def get_uniprots_for_compound(
    coconut_id, 
    threshold=0,
    filter_pa_pi=True,
    existing_conn=None):

    if isinstance(coconut_id, list) or isinstance(coconut_id, set):
        coconut_id = tuple(coconut_id)

    get_uniprots_sql = f'''
        SELECT c.coconut_id, t.target_name, u.acc
        FROM compounds AS c
        INNER JOIN activities AS a ON (c.compound_id=a.compound_id)
        INNER JOIN targets AS t ON (a.target_id=t.target_id)
        INNER JOIN targets_to_uniprot AS tu ON (a.target_id=tu.target_id)
        INNER JOIN uniprot AS u ON (tu.uniprot_id=u.uniprot_id)
        WHERE {f"c.coconut_id='{coconut_id}'" if isinstance(coconut_id, str)
            else f"c.coconut_id IN {coconut_id}"}
        {f"AND a.Pa>{threshold}" if threshold>0 else ""}
        {f"AND a.Pa>a.Pi" if filter_pa_pi else ""}
    '''

    return mysql_query(get_uniprots_sql, existing_conn=existing_conn)

def get_all_pathways(
    filter_actives=True,
    organisms=None, 
    existing_conn=None):
    if organisms is not None:
        if not isinstance(organisms, str):
            if isinstance(organisms, list) or isinstance(organisms, set):
                organisms = tuple(organisms)
            assert isinstance(organisms, tuple)
            print ("returning pathways for", len(organisms), "organisms")
            if len(organisms) == 1:
                organisms = organisms[0]
    active_filter = f'''
        INNER JOIN uniprot_to_pathway AS up ON (up.pathway_id=p.pathway_id)
        INNER JOIN targets_to_uniprot AS tu ON (up.uniprot_id=tu.uniprot_id)
        INNER JOIN targets AS t ON (t.target_id=tu.target_id)
    '''
    query = f'''
        SELECT DISTINCT p.pathway_name, p.organism
        FROM pathway AS p
        {active_filter if filter_actives else ""}
        {f"WHERE p.organism IN {organisms}" 
            if isinstance(organisms, tuple) else f"WHERE p.organism='{organisms}'" 
                if isinstance(organisms, str) else ""}
    '''
    return mysql_query(query, existing_conn=existing_conn)

def get_all_pathway_organisms(existing_conn=None):
    query = '''
        SELECT organism, COUNT(pathway_name)
        FROM pathway
        GROUP BY organism
    '''
    return mysql_query(query, existing_conn=existing_conn)

def get_all_reactions(
    filter_actives=True,
    organisms=None, 
    existing_conn=None):
    if organisms is not None:
        if not isinstance(organisms, str):
            if isinstance(organisms, list) or isinstance(organisms, set):
                organisms = tuple(organisms)
            assert isinstance(organisms, tuple)
            print ("returning pathways for", len(organisms), "organisms")
            if len(organisms) == 1:
                organisms = organisms[0]
    active_filter = f'''
        INNER JOIN uniprot_to_reaction AS ur ON (ur.reaction_id=r.reaction_id)
        INNER JOIN targets_to_uniprot AS tu ON (ur.uniprot_id=tu.uniprot_id)
        INNER JOIN targets AS t ON (t.target_id=tu.target_id)
    '''
    query = f'''
        SELECT DISTINCT r.reaction_name, r.organism
        FROM reaction AS r
        {active_filter if filter_actives else ""}
        {f"WHERE r.organism IN {organisms}" 
            if isinstance(organisms, tuple) else f"WHERE r.organism='{organisms}'" 
                if isinstance(organisms, str) else ""}
    '''
    return mysql_query(query, existing_conn=existing_conn)

def get_all_reaction_organisms(existing_conn=None):
    query = '''
        SELECT organism, COUNT(reaction_name)
        FROM reaction
        GROUP BY organism
    '''
    return mysql_query(query, existing_conn=existing_conn)

def get_all_pathways_for_compounds(
    coconut_ids,
    organism=None,
    filter_pa_pi=True,
    threshold=0,
    existing_conn=None,
    limit=None,
    ):
    assert filter_pa_pi
    if isinstance(coconut_ids, list) or isinstance(coconut_ids, set):
        coconut_ids = tuple(coconut_ids)
    else:
        assert isinstance(coconut_ids, str)
    # SELECT DISTINCT t.target_name, a.Pa, a.Pi, a.Pa-a.Pi, 
    # u.acc, tu.score, p.pathway_name, p.organism, up.evidence, p.pathway_url
    # INNER JOIN targets AS t ON (a.target_id=t.target_id)
    # {"AND a.Pa>a.Pi" if filter_pa_pi else ""}
    # {f"AND a.Pa>{threshold}" if threshold>0 else ""}

    query = f'''
        SELECT GROUP_CONCAT(DISTINCT(u.acc)) AS `Uniprot ACCS`, 
            COUNT(DISTINCT(u.acc)) AS `Number Uniprot ACCs`,
            p.pathway_name AS `Pathway Name`, p.organism AS `Organism`, 
            p.pathway_url AS `Pathway URL`
        FROM compounds AS c
        INNER JOIN activities AS a ON (c.compound_id=a.compound_id) 
        INNER JOIN targets_to_uniprot AS tu ON (a.target_id=tu.target_id)
        INNER JOIN uniprot AS u ON (tu.uniprot_id=u.uniprot_id)
        INNER JOIN uniprot_to_pathway AS up ON (tu.uniprot_id=up.uniprot_id)
        INNER JOIN pathway AS p ON (up.pathway_id=p.pathway_id)
        WHERE {f"c.coconut_id IN {coconut_ids}" if isinstance(coconut_ids, tuple)
            else f'c.coconut_id="{coconut_ids}"'}
        {f"AND a.above_{threshold}=(1)" if threshold>0 and filter_pa_pi else ""}
        {f"AND p.organism='{organism}'" if organism is not None else ""}
        GROUP BY p.pathway_name, p.organism
        {f"LIMIT {limit}" if limit is not None else ""}
    '''
    return mysql_query(query, existing_conn=existing_conn)

def get_all_reactions_for_compounds(
    coconut_ids,
    organism=None,
    filter_pa_pi=True,
    threshold=0,
    existing_conn=None,
    limit=None,
    ):
    assert filter_pa_pi
    if isinstance(coconut_ids, list):
        coconut_ids = tuple(coconut_ids)
    else:
        assert isinstance(coconut_ids, str)
    query = f'''
        SELECT GROUP_CONCAT(DISTINCT(u.acc)) AS `Uniprot ACCS`, 
            COUNT(DISTINCT(u.acc)) AS `Number Uniprot ACCs`,
            r.reaction_name AS `Reaction Name`, r.organism AS `Organism`, 
            r.reaction_url AS `Reaction URL`
        FROM compounds AS c
        INNER JOIN activities AS a ON (c.compound_id=a.compound_id) 
        INNER JOIN targets_to_uniprot AS tu ON (a.target_id=tu.target_id)
        INNER JOIN uniprot AS u ON (tu.uniprot_id=u.uniprot_id)
        INNER JOIN uniprot_to_reaction AS ur ON (tu.uniprot_id=ur.uniprot_id)
        INNER JOIN reaction AS r ON (ur.reaction_id=r.reaction_id)
        WHERE {f"c.coconut_id IN {coconut_ids}" if isinstance(coconut_ids, tuple)
            else f'c.coconut_id="{coconut_ids}"'}
        {f"AND a.above_{threshold}=(1)" if threshold>0 and filter_pa_pi else ""}
        {f"AND r.organism='{organism}'" if organism is not None else ""}
        GROUP BY r.reaction_name, r.organism
        {f"LIMIT {limit}" if limit is not None else ""}
    '''
    return mysql_query(query, existing_conn=existing_conn)


if __name__ == "__main__":

    records = get_all_reactions(organisms="Homo sapiens")

    print (len(records))
    for record in records[:5]:
        print (record)