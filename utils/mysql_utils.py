import sys
import os.path
sys.path.insert(1, 
    os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))

import mysql.connector
from mysql.connector.errors import ProgrammingError, IntegrityError

from utils.io import load_json
from utils.pass_utils import remove_invalid_characters, parse_pass_spectra, get_all_targets, get_categories, get_targets_for_category, get_all_compounds

from functools import partial

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

def create_tables(): 

    # compounds = get_all_compounds()
    # targets = get_all_targets()
    # # n_targets = len(targets)

    # from natural_products.backend import get_multiple_compound_info
    # # compound_info = get_multiple_compound_info()

    # conn = connect_to_mysqldb()

    # # create compounds table
    # create_compounds_table = f"CREATE TABLE `compounds` (compound_id MEDIUMINT, coconut_id VARCHAR(255), name VARCHAR(2000), "+\
    #     "formula VARCHAR(255), smiles VARCHAR(1000), clean_smiles VARCHAR(1000), PRIMARY KEY (compound_id))"
    # mysql_create_table(create_compounds_table, existing_conn=conn)

    # # insert into compounds table
    # insert_compounds_sql = "INSERT INTO `compounds`" + " (compound_id, coconut_id, name, formula, smiles) VALUES (%s, %s, %s, %s, %s) "+\
    #     "ON DUPLICATE KEY UPDATE compound_id=VALUES(compound_id), coconut_id=VALUES(coconut_id), name=VALUES(name), formula=VALUES(formula), smiles=VALUES(smiles)"
    # rows = [(i, c, *compound_info[c]) for c, i in compounds.items()]
    # mysql_insert_many(insert_compounds_sql, rows, existing_conn=conn)

    # # create targets table
    # create_targets_table = f"CREATE TABLE `targets` (target_id SMALLINT, target_name VARCHAR(255), PRIMARY KEY (target_id))"
    # mysql_create_table(create_targets_table, existing_conn=conn)

    # # insert into targets table
    # insert_targets_sql = "INSERT INTO `targets`" + " (target_id, target_name) VALUES (%s, %s) ON DUPLICATE KEY UPDATE target_id=target_id"
    # rows = [(i, t) for t, i in targets.items()]
    # mysql_insert_many(insert_targets_sql, rows, existing_conn=conn)

    # # # create activities table
    # create_activities_table = f"CREATE TABLE `activities` (compound_id MEDIUMINT, target_id SMALLINT, Pa SMALLINT, Pi SMALLINT, " +\
    #     "PRIMARY KEY (compound_id, target_id), FOREIGN KEY (compound_id) REFERENCES compounds(compound_id), FOREIGN KEY (target_id) REFERENCES targets(target_id))"
    # mysql_create_table(create_activities_table, existing_conn=conn)

    # # create categories table
    # create_categories_table = f"CREATE TABLE `categories` (category_id TINYINT, category_name VARCHAR(255), " +\
    #     "PRIMARY KEY (category_id))"
    # mysql_create_table(create_categories_table, existing_conn=conn)

    # # insert into categories table
    # categories = {c: i for i, c in enumerate(sorted(get_categories()))}
    # insert_categories_sql = "INSERT INTO `categories`" + " (category_id, category_name) VALUES (%s, %s) ON DUPLICATE KEY UPDATE category_id=category_id"
    # rows = [(i, c) for c, i in categories.items()]
    # mysql_insert_many(insert_categories_sql, rows, existing_conn=conn)

    # # create categories members table
    # create_categories_members_table = f"CREATE TABLE `category_members` (category_id TINYINT, target_id SMALLINT, " +\
    #     "PRIMARY KEY (category_id, target_id), "+\
    #         "FOREIGN KEY (category_id) REFERENCES categories(category_id), "+\
    #         "FOREIGN KEY (target_id) REFERENCES targets(target_id) )"
    # mysql_create_table(create_categories_members_table, existing_conn=conn)

    # # insert into categories table
    # insert_categories_members_sql = "INSERT INTO `category_members`" + " (category_id, target_id) VALUES (%s, %s) "+\
    #     "ON DUPLICATE KEY UPDATE category_id=category_id, target_id=target_id"

    # for c, i in categories.items():
    #     category_targets = get_targets_for_category(c)
    #     rows = [(i, t_i) for t, t_i in category_targets.items()]
    #     mysql_insert_many(insert_categories_members_sql, rows, existing_conn=conn)

    # create_uniprot_table = '''
    #     CREATE TABLE uniprot (
    #         uniprot_id MEDIUMINT NOT NULL AUTO_INCREMENT,
    #         acc VARCHAR(10) NOT NULL UNIQUE,
    #         PRIMARY KEY(id)
    #     )
    # '''

    # mysql_create_table(create_uniprot_table)

    create_target_uniprot_link_table = '''
    CREATE TABLE targets_to_uniprot(
        target_id SMALLINT, 
        uniprot_id MEDIUMINT,
        score FLOAT,
        PRIMARY KEY(target_id, uniprot_id),
        FOREIGN KEY(target_id) REFERENCES targets(target_id),
        FOREIGN KEY(uniprot_id) REFERENCES uniprot(uniprot_id)
    )
    '''

    mysql_create_table(create_target_uniprot_link_table)

    # create_pathway_table = '''
    #     CREATE TABLE pathway (
    #         pathway_id SMALLINT NOT NULL AUTO_INCREMENT,
    #         reactome_identifier VARCHAR(255) NOT NULL UNIQUE,
    #         pathway_name VARCHAR(255) NOT NULL,
    #         organism VARCHAR(255) NOT NULL,
    #         pathway_url VARCHAR(255),
    #         PRIMARY KEY(pathway_id)
    #     )
    # '''

    # mysql_create_table(create_pathway_table)

    # create_uniprot_to_pathway_table = '''
    #     CREATE TABLE uniprot_to_pathway (
    #         uniprot_id MEDIUMINT NOT NULL,
    #         pathway_id SMALLINT NOT NULL,
    #         evidence VARCHAR(25),
    #         PRIMARY KEY(uniprot_id, pathway_id),
    #         FOREIGN KEY(uniprot_id) REFERENCES uniprot(uniprot_id),
    #         FOREIGN KEY(pathway_id) REFERENCES pathway(pathway_id)
    #     )
    # '''

    # mysql_create_table(create_uniprot_to_pathway_table)

    # create_reaction_table = '''
    #     CREATE TABLE reaction (
    #         reaction_id MEDIUMINT NOT NULL AUTO_INCREMENT,
    #         reactome_identifier VARCHAR(255) NOT NULL UNIQUE,
    #         reaction_name VARCHAR(255) NOT NULL,
    #         organism VARCHAR(255) NOT NULL,
    #         reaction_url VARCHAR(255),
    #         PRIMARY KEY(reaction_id)
    #     )
    # '''

    # mysql_create_table(create_reaction_table)

    # create_uniprot_to_reaction_table = '''
    #     CREATE TABLE uniprot_to_reaction (
    #         uniprot_id MEDIUMINT NOT NULL,
    #         reaction_id MEDIUMINT NOT NULL,
    #         evidence VARCHAR(25),
    #         PRIMARY KEY(uniprot_id, reaction_id),
    #         FOREIGN KEY(uniprot_id) REFERENCES uniprot(uniprot_id),
    #         FOREIGN KEY(reaction_id) REFERENCES reaction(reaction_id)
    #     )
    # '''

    # mysql_create_table(create_uniprot_to_reaction_table)

    return 0


def get_all_targets_and_categories(categories=None):

    if categories is not None:
        if isinstance(categories, list) or isinstance(categories, set):
            categories = tuple(categories)

    query = f'''
       SELECT c.category_name, t.target_name
       FROM categories AS c
       INNER JOIN category_members AS m ON (c.category_id=m.category_id) 
       INNER JOIN targets AS t ON (m.target_id=t.target_id)
       {f"WHERE c.category_name IN {categories}" if categories is not None else ""}
    '''

    return mysql_query(query)


def migrate_SDF_to_mysql(sdf_file):
    '''
        read in PASS predictions from SDF file and add to SQL database
        only valid for GUI PASS version
    '''
    print ("reading SDF file from", sdf_file, "and inserting into SQL database")

    compound_to_id = get_all_compounds()
    target_to_id = get_all_targets()

    categories = get_categories()

    import pandas as pd

    # from rdkit.Chem.PandasTools import LoadSDF

    from collections import defaultdict

    from utils.rdkit_utils import LoadSDF

    chunks = LoadSDF(sdf_file, smilesName="SMILES", molColName=None, chunksize=5000)

    sql = '''
        INSERT INTO activities (compound_id, target_id, Pa, Pi) 
            VALUES (%s, %s, %s, %s) 
        ON DUPLICATE KEY UPDATE compound_id=compound_id, target_id=target_id
        '''

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

            rows = (
                (compound, target_id, *row) 
                for target_id in targets
                for compound, row in zip(compound_ids, targets[target_id])
            )

            mysql_insert_many(sql, rows, )

            print ("completed category", category, "for chunk no", chunk_no,
                "for file", sdf_file)
            print ("################")
            print ()
       
        print ("completed chunk number", chunk_no, "for file", sdf_file)
        print ("################")
        print ()


def add_clean_smiles():

    from natural_products.backend import get_multiple_compound_info

    smiles = get_multiple_compound_info(
        columns=("compound_id", "smiles"))

    from rdkit import Chem
    from utils.rdkit_utils import standardise_smi

    standardised_smiles = list(
        map(lambda smi: standardise_smi(smi, return_smiles=True), 
            (smi for c, smi in smiles)))

    # assert len(smiles) == len(standardised_smiles)

    insert_standard_smiles_sql = '''
        INSERT INTO compounds (compound_id, clean_smiles)
            VALUES (%s, %s)
            ON DUPLICATE KEY UPDATE compound_id=VALUES(compound_id), clean_smiles=VALUES(clean_smiles)
        '''
    rows = (
        (s, standardise_smi)
            for (s, smi), standardise_smi in zip(smiles, standardised_smiles)
    )
    mysql_insert_many(insert_standard_smiles_sql, rows)

def add_uniprot_accs(uniprot_accs, existing_conn=None):
    print ("inserting", len(uniprot_accs), "UNIPROT ACCs")
    insert_uniprots_sql = '''
        INSERT INTO uniprot (acc) VALUES (%s)
            ON DUPLICATE KEY UPDATE acc=acc
    '''
    mysql_insert_many(insert_uniprots_sql, 
        ((uniprot_acc,) for uniprot_acc in uniprot_accs),
        existing_conn=existing_conn)

    return 0

def add_target_to_uniprot(chunksize=50):
    from utils.genenames_utils import search_for_targets

    '''
        MECHANISMS is (too?) general
        
        COMPLETED
        EFFECTS
        ANTITARGETS
        GENE EXPRESSIONS
        TRANSPORTERS
        TOXICITY
        METABOLISM
    '''

    # get valid targets
    get_target_names_sql = f'''
        SELECT t.target_id, t.target_name
        FROM targets AS t
        INNER JOIN category_members AS cm ON (t.target_id=cm.target_id)
        INNER JOIN categories AS c ON (cm.category_id=c.category_id)
        WHERE t.target_id NOT IN (
            SELECT DISTINCT target_id FROM targets_to_uniprot)
        AND c.category_name="METABOLISM"
    '''
    targets = mysql_query(get_target_names_sql)
    target_map = {t:i for i, t in targets}

    n_targets = len(targets)
    n_chunks = n_targets // chunksize + 1
    assert n_chunks * chunksize >= n_targets

    for chunk_no in range(n_chunks):

        targets_to_uniprots = search_for_targets(
            [target_name for target_id, target_name in targets[chunk_no*chunksize:(chunk_no+1)*chunksize]],
            key="uniprot_ids")

        unique_uniprots = {uniprot
            for target in targets_to_uniprots
            for uniprot, score in targets_to_uniprots[target]}

        add_uniprot_accs(unique_uniprots)

        # get all uniprots for mapping (reduce queries)
        get_all_uniprots_sql = '''
            SELECT uniprot_id, acc
            from uniprot
        '''
        uniprot_map = mysql_query(get_all_uniprots_sql)
        uniprot_map = {u: i for i, u in uniprot_map}

        insert_targets_to_uniprot_sql = '''
            INSERT INTO targets_to_uniprot (target_id, uniprot_id, score)
                VALUES (%s, %s, %s) ON DUPLICATE KEY UPDATE
                target_id=target_id, uniprot_id=uniprot_id, score=score
        '''

        rows = []

        for i, (target_name, uniprots) in \
            enumerate(targets_to_uniprots.items()):

            target_rows = [
                (target_map[target_name], uniprot_map[uniprot], score) # target and uniprot must be in database
                    for uniprot, score in uniprots
            ]
            assert len(target_rows) == len(uniprots)
            rows.extend(target_rows)

            print ("completed target", i+1, "/", len(targets_to_uniprots))

        mysql_insert_many(insert_targets_to_uniprot_sql, rows)
        print ("completed chunk no", chunk_no+1, "/", n_chunks)

def get_uniprots_for_targets(targets, existing_conn=None):
    if not isinstance(targets, tuple):
        assert isinstance(targets, list) or isinstance(targets, set)
        targets = tuple(targets)
    print ("querying mySQL database for UNIPROTS associated with",
        len(targets), "targets")
    query = f'''
        SELECT DISTINCT t.target_name, u.acc, tu.score
        FROM targets AS t
        INNER JOIN targets_to_uniprot AS tu ON (t.target_id=tu.target_id)
        INNER JOIN uniprot AS u ON (tu.uniprot_id=u.uniprot_id)
        WHERE t.target_name IN {targets}
    '''
    return mysql_query(query, existing_conn=existing_conn)

def get_uniprots_for_compound(
    coconut_id, 
    threshold=0,
    filter_pa_pi=True):

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

    return mysql_query(get_uniprots_sql)

def add_pathways_and_reactions():
    import pandas as pd

    # read pathways
    pathways = pd.read_csv("/home/david/Desktop/all_pathways.csv",
        index_col=0)

    insert_pathways_sql = '''
        INSERT INTO pathway (reactome_identifier,
            pathway_name, organism, pathway_url)
            VALUES (%s, %s, %s, %s)
            ON DUPLICATE KEY UPDATE 
            reactome_identifier=reactome_identifier
    '''

    rows = (
        (row["identifier"], row["pathway"],
            row["organism"], row["pathway_url"])
            for _, row in pathways.iterrows() 
    )

    mysql_insert_many(insert_pathways_sql, rows)

    # read reactions
    reactions = pd.read_csv("/home/david/Desktop/reactions.csv",
        index_col=0)

    insert_reactions_sql = '''
        INSERT INTO reaction (reactome_identifier,
            reaction_name, organism, reaction_url)
            VALUES (%s, %s, %s, %s)
            ON DUPLICATE KEY UPDATE 
            reactome_identifier=reactome_identifier
    '''

    rows = (
        (row["identifier"], row["reaction"],
            row["organism"], row["reaction_url"])
            for _, row in reactions.iterrows() 
    )

    mysql_insert_many(insert_reactions_sql, rows)

def add_uniprot_to_pathways():
    import pandas as pd

    uniprot_to_pathways = pd.read_csv(
        "/home/david/Desktop/uniprot_to_pathway.csv", 
        index_col=0)

    conn = connect_to_mysqldb()

    unique_uniprots = set(uniprot_to_pathways["acc"]) # add all associated uniprots to database
    add_uniprot_accs(unique_uniprots, existing_conn=conn)

    uniprot_query = '''
        SELECT uniprot_id, acc
        FROM uniprot
    '''
    records = mysql_query(uniprot_query, existing_conn=conn)
    uniprot_to_id = {
        a: i for i, a in records
    }

    pathway_query = '''
        SELECT pathway_id, pathway_name, organism
        FROM pathway
    '''
    records = mysql_query(pathway_query, existing_conn=conn)
    pathway_to_id = {
        (n, o): i for i, n, o in records
    }

    rows = []
    for _, row in uniprot_to_pathways.iterrows():
        accession_number = row["acc"]
        pathway_name = row["pathway"]
        organism = row["organism"]
        evidence = row["evidence"]

        rows.append((uniprot_to_id[accession_number],
            pathway_to_id[(pathway_name, organism)],
            evidence))

    insert_uniprot_to_pathways_sql = '''
        INSERT INTO uniprot_to_pathway 
            (uniprot_id, pathway_id, evidence)
            VALUES(%s, %s, %s)
            ON DUPLICATE KEY UPDATE uniprot_id=uniprot_id
    '''

    mysql_insert_many(insert_uniprot_to_pathways_sql, rows,
        existing_conn=conn)

def add_uniprot_to_reactions():
    import pandas as pd

    uniprot_to_reaction = pd.read_csv(
        "/home/david/Desktop/uniprot_to_reaction.csv", 
        index_col=0)

    conn = connect_to_mysqldb()

    unique_uniprots = set(uniprot_to_reaction["acc"])
    add_uniprot_accs(unique_uniprots, existing_conn=conn)

    uniprot_query = '''
        SELECT uniprot_id, acc
        FROM uniprot
    '''
    records = mysql_query(uniprot_query, existing_conn=conn)
    uniprot_to_id = {
        a: i for i, a in records
    }

    reaction_query = '''
        SELECT reaction_id, reaction_name, organism
        FROM reaction
    '''
    records = mysql_query(reaction_query, existing_conn=conn)
    reaction_to_id = {
        (n, o): i for i, n, o in records
    }

    rows = []
    for _, row in uniprot_to_reaction.iterrows():
        accession_number = row["acc"]
        reaction_name = row["reaction"]
        organism = row["organism"]
        evidence = row["evidence"]

        rows.append((uniprot_to_id[accession_number],
            reaction_to_id[(reaction_name, organism)],
            evidence))

    insert_uniprot_to_reaction_sql = '''
        INSERT INTO uniprot_to_reaction
            (uniprot_id, reaction_id, evidence)
            VALUES(%s, %s, %s)
            ON DUPLICATE KEY UPDATE uniprot_id=uniprot_id
    '''

    mysql_insert_many(insert_uniprot_to_reaction_sql, rows,
        existing_conn=conn)

def get_all_pathways(
    actives=True,
    organism=None, 
    existing_conn=None):
    query = f'''
        SELECT pathway_name, organism
        FROM pathway
        {"WHERE pathway_id IN (SELECT DISTINCT pathway_id from uniprot_to_pathway)"
            if actives else ""}
    '''
    if organism is not None:
        query += f"AND organism=\"{organism}\""
    return mysql_query(query, existing_conn=existing_conn)

def get_all_pathways_for_compounds(
    coconut_ids,
    organism=None,
    filter_pa_pi=True,
    threshold=0,
    existing_conn=None,
    limit=None,
    ):
    if isinstance(coconut_ids, list) or isinstance(coconut_ids, set):
        coconut_ids = tuple(coconut_ids)
    else:
        assert isinstance(coconut_ids, str)
    filter_pa_pi = False
    query = f'''
        SELECT DISTINCT t.target_name, a.Pa, a.Pi, a.Pa-a.Pi, 
            u.acc, tu.score, p.pathway_name, p.organism, up.evidence, p.pathway_url
        FROM compounds AS c
        INNER JOIN activities AS a ON (c.compound_id=a.compound_id) 
        INNER JOIN targets AS t ON (a.target_id=t.target_id)
        INNER JOIN targets_to_uniprot AS tu ON (a.target_id=tu.target_id)
        INNER JOIN uniprot AS u ON (tu.uniprot_id=u.uniprot_id)
        INNER JOIN uniprot_to_pathway AS up ON (tu.uniprot_id=up.uniprot_id)
        INNER JOIN pathway AS p ON (up.pathway_id=p.pathway_id)
        WHERE {f"c.coconut_id IN {coconut_ids}" if isinstance(coconut_ids, tuple)
            else f'c.coconut_id="{coconut_ids}"'}
        {"AND a.Pa>a.Pi" if filter_pa_pi else ""}
        {f"AND a.Pa>{threshold}" if threshold>0 else ""}
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
    if isinstance(coconut_ids, list):
        coconut_ids = tuple(coconut_ids)
    else:
        assert isinstance(coconut_ids, str)
    filter_pa_pi = False
    query = f'''
        SELECT DISTINCT t.target_name, a.Pa, a.Pi, a.Pa-a.Pi, 
            u.acc, tu.score, r.reaction_name, r.organism, ur.evidence, r.reaction_url
        FROM compounds AS c
        INNER JOIN activities AS a ON (c.compound_id=a.compound_id) 
        INNER JOIN targets AS t ON (a.target_id=t.target_id) 
        INNER JOIN targets_to_uniprot AS tu ON (a.target_id=tu.target_id)
        INNER JOIN uniprot AS u ON (tu.uniprot_id=u.uniprot_id)
        INNER JOIN uniprot_to_reaction AS ur ON (tu.uniprot_id=ur.uniprot_id)
        INNER JOIN reaction AS r ON (ur.reaction_id=r.reaction_id)
        WHERE {f"c.coconut_id IN {coconut_ids}" if isinstance(coconut_ids, tuple)
            else f'c.coconut_id="{coconut_ids}"'}
        {"AND a.Pa>a.Pi" if filter_pa_pi else ""}
        {f"AND a.Pa>{threshold}" if threshold>0 else ""}
        {f"LIMIT {limit}" if limit is not None else ""}
    '''
    return mysql_query(query, existing_conn=existing_conn)


def get_all_reactions(
    actives=True,
    organism=None, 
    existing_conn=None):
    query = f'''
        SELECT reaction_name, organism
        FROM reaction
        {"WHERE reaction_id IN (SELECT DISTINCT reaction_id from uniprot_to_reaction)"
            if actives else ""}
    '''
    if organism is not None:
        query += f"AND organism=\"{organism}\""
    return mysql_query(query, existing_conn=existing_conn)


if __name__ == "__main__":
    # create_tables()
    # add_target_to_uniprot()
    # pass

    # records = get_all_reactions_for_compounds(coconut_ids="CNP0000002", organism="Homo sapiens", threshold=950)

    # for record in records:
    #     print (record)

    # records = get_uniprots_for_targets(targets=("Diarrhea", "Yawning"))
    records = get_all_reactions_for_compounds("CNP0000002", threshold=900)
    for record in records[:10]:
        print (record)

    # pathway_name = "Urea cycle"
    # reaction_name = "Expression of ABCA1"
    # organism = "Homo sapiens"

    # records = get_compounds_for_reaction(reaction_name, organism,
    #     threshold=0, filter_pa_pi=True, include_targets=True)

    # print (len(records))
    # # print (len(set((record[0] for record in records))))
    # for record in records[:10]:
    #     print (record)