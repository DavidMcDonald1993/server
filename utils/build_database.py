import sys
import os.path
sys.path.insert(1, 
    os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))

from utils.mysql_utils import mysql_create_table, mysql_insert_many, mysql_query, connect_to_mysqldb
from utils.pass_utils import (remove_invalid_characters, parse_pass_spectra, get_all_targets, 
    get_categories, get_targets_for_category, get_all_compounds)

from natural_products.backend import get_coconut_compound_info_from_mongo
from utils.queries import get_info_for_multiple_compounds

from functools import partial

def create_and_populate_compounds_table(): 

    compounds = get_all_compounds()
    compound_info = get_coconut_compound_info_from_mongo()

    conn = connect_to_mysqldb()

    # create compounds table
    create_compounds_table_sql = f'''
    CREATE TABLE `compounds` (
        `compound_id` mediumint NOT NULL,
        `coconut_id` varchar(15) DEFAULT NULL,
        `name` varchar(1200) DEFAULT NULL,
        `formula` varchar(50) DEFAULT NULL,
        `smiles` varchar(600) DEFAULT NULL,
        `clean_smiles` varchar(600) DEFAULT NULL,
        `image_path` varchar(50) GENERATED ALWAYS AS (concat(_utf8mb4'compound_images',_utf8mb4'/',cast((`compound_id` DIV 1024) as char charset utf8mb4),_utf8mb4'/',`compound_id`,_utf8mb4'.png')) VIRTUAL NOT NULL,
        PRIMARY KEY (`compound_id`),
        KEY `coconut_id` (`coconut_id`)
        ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci
    '''
    mysql_create_table(create_compounds_table_sql, existing_conn=conn)

    # insert into compounds table
    insert_compounds_sql = '''
    INSERT INTO `compounds` (compound_id, coconut_id, name, formula, smiles) 
        VALUES (%s, %s, %s, %s, %s)
        ON DUPLICATE KEY UPDATE compound_id=VALUES(compound_id), 
            coconut_id=VALUES(coconut_id), name=VALUES(name), formula=VALUES(formula), smiles=VALUES(smiles)
    '''
    rows = [(i, c, *compound_info[c]) for c, i in compounds.items()]
    mysql_insert_many(insert_compounds_sql, rows, existing_conn=conn)

    return 0

def add_clean_smiles_to_compounds():

    smiles = get_info_for_multiple_compounds(
        columns=("compound_id", "smiles"))

    from utils.rdkit_utils import standardise_smi

    standardised_smiles = list(
        map(lambda smi: standardise_smi(smi, return_smiles=True), 
            (smi for c, smi in smiles)))

    insert_standard_smiles_sql = '''
    INSERT INTO compounds (compound_id, clean_smiles)
        VALUES (%s, %s)
        ON DUPLICATE KEY UPDATE compound_id=VALUES(compound_id), clean_smiles=VALUES(clean_smiles)
    '''
    rows = (
        (c, standardised_smi)
            for (c, smi), standardised_smi in zip(smiles, standardised_smiles)
    )
    return mysql_insert_many(insert_standard_smiles_sql, rows)

def create_and_populate_targets_table():

    conn = connect_to_mysqldb()

    # create targets table
    create_targets_table_sql = '''
    CREATE TABLE `targets` (
        `target_id` smallint NOT NULL,
        `target_name` varchar(115) DEFAULT NULL,
        PRIMARY KEY (`target_id`),
        KEY `target_name` (`target_name`)
        ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci
    '''
    mysql_create_table(create_targets_table_sql, existing_conn=conn)

    targets = get_all_targets()

    # insert into targets table
    insert_targets_sql = '''
    INSERT INTO `targets` (target_id, target_name) 
        VALUES (%s, %s) ON DUPLICATE KEY UPDATE target_id=target_id
    '''
    rows = [(i, t) for t, i in targets.items()]
    mysql_insert_many(insert_targets_sql, rows, existing_conn=conn)

    return 0

def create_activities_table():
    create_activities_table_sql = '''
    CREATE TABLE `activities` (
    `compound_id` mediumint NOT NULL,
    `target_id` smallint NOT NULL,
    `Pa` smallint DEFAULT NULL,
    `Pi` smallint DEFAULT NULL,
    `above_950` bit(1) GENERATED ALWAYS AS (((`Pa` > `Pi`) and (`Pa` > 950))) VIRTUAL NOT NULL,
    `above_900` bit(1) GENERATED ALWAYS AS (((`Pa` > `Pi`) and (`Pa` > 900))) VIRTUAL NOT NULL,
    `above_850` bit(1) GENERATED ALWAYS AS (((`Pa` > `Pi`) and (`Pa` > 850))) VIRTUAL NOT NULL,
    `above_800` bit(1) GENERATED ALWAYS AS (((`Pa` > `Pi`) and (`Pa` > 800))) VIRTUAL NOT NULL,
    `above_750` bit(1) GENERATED ALWAYS AS (((`Pa` > `Pi`) and (`Pa` > 750))) VIRTUAL NOT NULL,
    `above_700` bit(1) GENERATED ALWAYS AS (((`Pa` > `Pi`) and (`Pa` > 700))) VIRTUAL NOT NULL,
    `above_650` bit(1) GENERATED ALWAYS AS (((`Pa` > `Pi`) and (`Pa` > 650))) VIRTUAL NOT NULL,
    `above_600` bit(1) GENERATED ALWAYS AS (((`Pa` > `Pi`) and (`Pa` > 600))) VIRTUAL NOT NULL,
    `above_550` bit(1) GENERATED ALWAYS AS (((`Pa` > `Pi`) and (`Pa` > 550))) VIRTUAL NOT NULL,
    `above_500` bit(1) GENERATED ALWAYS AS (((`Pa` > `Pi`) and (`Pa` > 500))) VIRTUAL NOT NULL,
    `above_450` bit(1) GENERATED ALWAYS AS (((`Pa` > `Pi`) and (`Pa` > 450))) VIRTUAL NOT NULL,
    `above_400` bit(1) GENERATED ALWAYS AS (((`Pa` > `Pi`) and (`Pa` > 400))) VIRTUAL NOT NULL,
    `above_350` bit(1) GENERATED ALWAYS AS (((`Pa` > `Pi`) and (`Pa` > 350))) VIRTUAL NOT NULL,
    `above_300` bit(1) GENERATED ALWAYS AS (((`Pa` > `Pi`) and (`Pa` > 300))) VIRTUAL NOT NULL,
    `above_250` bit(1) GENERATED ALWAYS AS (((`Pa` > `Pi`) and (`Pa` > 250))) VIRTUAL NOT NULL,
    `above_200` bit(1) GENERATED ALWAYS AS (((`Pa` > `Pi`) and (`Pa` > 200))) VIRTUAL NOT NULL,
    `above_150` bit(1) GENERATED ALWAYS AS (((`Pa` > `Pi`) and (`Pa` > 150))) VIRTUAL NOT NULL,
    `above_100` bit(1) GENERATED ALWAYS AS (((`Pa` > `Pi`) and (`Pa` > 100))) VIRTUAL NOT NULL,
    `above_50` bit(1) GENERATED ALWAYS AS (((`Pa` > `Pi`) and (`Pa` > 50))) VIRTUAL NOT NULL,
    `confidence_score` smallint GENERATED ALWAYS AS ((`Pa` - `Pi`)) VIRTUAL NOT NULL,
    PRIMARY KEY (`compound_id`,`target_id`),
    KEY `target_id` (`target_id`),
    KEY `above_950_target_id` (`above_950`,`target_id`),
    KEY `above_900_target_id` (`above_900`,`target_id`),
    KEY `above_850_target_id` (`above_850`,`target_id`),
    KEY `above_800_target_id` (`above_800`,`target_id`),
    KEY `above_750_target_id` (`above_750`,`target_id`),
    KEY `above_700_target_id` (`above_700`,`target_id`),
    KEY `above_650_target_id` (`above_650`,`target_id`),
    CONSTRAINT `activities_ibfk_1` FOREIGN KEY (`compound_id`) REFERENCES `compounds` (`compound_id`),
    CONSTRAINT `activities_ibfk_2` FOREIGN KEY (`target_id`) REFERENCES `targets` (`target_id`)
    ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci
    '''
    mysql_create_table(create_activities_table_sql, )

    return 0

def populate_activities_from_PASS_sdf_file(sdf_file):
    '''
        read in PASS predictions from SDF file and add to SQL database
        only valid for GUI PASS version
    '''
    print ("reading SDF file from", sdf_file, "and inserting into SQL database")

    compound_to_id = get_all_compounds()
    target_to_id = get_all_targets()
    categories = get_categories()

    import pandas as pd
    from collections import defaultdict
    # from rdkit.Chem.PandasTools import LoadSDF
    from utils.rdkit_utils import LoadSDF

    chunks = LoadSDF(sdf_file, smilesName="SMILES", molColName=None, chunksize=5000)

    insert_activities_sql = '''
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

            mysql_insert_many(insert_activities_sql, rows, )

            print ("completed category", category, "for chunk no", chunk_no,
                "for file", sdf_file)
            print ("################")
            print ()
       
        print ("completed chunk number", chunk_no, "for file", sdf_file)
        print ("################")
        print ()


def create_and_populate_categories_table():

    conn = connect_to_mysqldb()

    # create categories table
    create_categories_table_sql = '''
    CREATE TABLE `categories` (
    `category_id` tinyint NOT NULL,
    `category_name` varchar(20) DEFAULT NULL,
    PRIMARY KEY (`category_id`)
    ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci
    '''
    mysql_create_table(create_categories_table_sql, existing_conn=conn)

    # insert into categories table
    categories = {c: i 
        for i, c in enumerate(sorted(get_categories()))}
    insert_categories_sql = '''
    INSERT INTO `categories` (category_id, category_name) 
        VALUES (%s, %s) 
        ON DUPLICATE KEY UPDATE category_id=category_id
    '''
    rows = [(i, c) for c, i in categories.items()]
    mysql_insert_many(insert_categories_sql, rows, existing_conn=conn)

    return 0

def create_and_populate_category_members():

    conn = connect_to_mysqldb()

    # create categories members table
    create_categories_members_table_sql = '''
    CREATE TABLE `category_members` (
        category_id TINYINT, 
        target_id SMALLINT,
        PRIMARY KEY (category_id, target_id),
        FOREIGN KEY (category_id) REFERENCES categories(category_id),
        FOREIGN KEY (target_id) REFERENCES targets(target_id) 
    )
    '''
    mysql_create_table(create_categories_members_table_sql, existing_conn=conn)

    # insert into categories members table
    insert_categories_members_sql = '''
        INSERT INTO `category_members` (category_id, target_id) 
        VALUES (%s, %s)
        ON DUPLICATE KEY UPDATE category_id=category_id, target_id=target_id
    '''
    for i, c in enumerate(sorted(get_categories())):
        category_targets = get_targets_for_category(c)
        rows = [(i, t_i) for t, t_i in category_targets.items()]
        mysql_insert_many(insert_categories_members_sql, rows, existing_conn=conn)

    return 0


def create_uniprot_table():
    create_uniprot_table_sql = '''
    CREATE TABLE uniprot (
        uniprot_id MEDIUMINT NOT NULL AUTO_INCREMENT,
        acc VARCHAR(10) NOT NULL UNIQUE,
        PRIMARY KEY(uniprot_id),
        UNIQUE KEY `acc` (`acc`)
    )
    '''

    return mysql_create_table(create_uniprot_table_sql)

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

def create_target_to_uniprot_table():

    create_target_to_uniprot_table_sql = '''
    CREATE TABLE targets_to_uniprot(
        target_id SMALLINT, 
        uniprot_id MEDIUMINT,
        score FLOAT,
        PRIMARY KEY(target_id, uniprot_id),
        FOREIGN KEY(target_id) REFERENCES targets(target_id),
        FOREIGN KEY(uniprot_id) REFERENCES uniprot(uniprot_id)
    )
    '''

    return mysql_create_table(create_target_to_uniprot_table_sql)

def populate_targets_to_uniprot(chunksize=50):
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
            SELECT DISTINCT target_id 
            FROM targets_to_uniprot
        )
    '''
    targets = mysql_query(get_target_names_sql)
    target_map = {t: i for i, t in targets}

    n_targets = len(targets)
    n_chunks = n_targets // chunksize + 1
    assert n_chunks * chunksize >= n_targets

    for chunk_no in range(n_chunks):

        targets_to_uniprots = search_for_targets(
            search_terms=[target_name 
                for target_id, target_name in targets[chunk_no*chunksize:(chunk_no+1)*chunksize]],
            key="uniprot_ids")

        unique_uniprots = {uniprot
            for target in targets_to_uniprots
            for uniprot, score in targets_to_uniprots[target]
        }

        add_uniprot_accs(unique_uniprots) # add accs not in database

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

            print ("completed target", i+1, "/", len(targets_to_uniprots),
                "for chunk number", chunk_no+1, "/", n_chunks)

        mysql_insert_many(insert_targets_to_uniprot_sql, rows)
        print ("completed chunk no", chunk_no+1, "/", n_chunks)
        print ("#"*50)
        print ()

def create_pathway_table():

    create_pathway_table_sql = '''
        CREATE TABLE pathway (
            pathway_id SMALLINT NOT NULL AUTO_INCREMENT,
            reactome_identifier VARCHAR(255) NOT NULL UNIQUE,
            pathway_name VARCHAR(255) NOT NULL,
            organism VARCHAR(255) NOT NULL,
            pathway_url VARCHAR(255),
            PRIMARY KEY(pathway_id)
        )
    '''

    return mysql_create_table(create_pathway_table_sql)


def populate_pathways():
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

    return mysql_insert_many(insert_pathways_sql, rows)

def populate_reactions():
    import pandas as pd

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

    return mysql_insert_many(insert_reactions_sql, rows)

def create_uniprot_to_pathway_table():

    create_uniprot_to_pathway_table_sql = '''
        CREATE TABLE uniprot_to_pathway (
            uniprot_id MEDIUMINT NOT NULL,
            pathway_id SMALLINT NOT NULL,
            evidence VARCHAR(25),
            PRIMARY KEY(uniprot_id, pathway_id),
            FOREIGN KEY(uniprot_id) REFERENCES uniprot(uniprot_id),
            FOREIGN KEY(pathway_id) REFERENCES pathway(pathway_id)
        )
    '''

    return mysql_create_table(create_uniprot_to_pathway_table_sql)

def populate_uniprot_to_pathways():
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

        rows.append(
            (uniprot_to_id[accession_number], 
                pathway_to_id[(pathway_name, organism)], evidence)
        )

    insert_uniprot_to_pathways_sql = '''
        INSERT INTO uniprot_to_pathway 
            (uniprot_id, pathway_id, evidence)
            VALUES(%s, %s, %s)
            ON DUPLICATE KEY UPDATE uniprot_id=uniprot_id
    '''

    return mysql_insert_many(insert_uniprot_to_pathways_sql, rows, existing_conn=conn)

def create_reaction_table():

    create_reaction_table_sql = '''
        CREATE TABLE reaction (
            reaction_id MEDIUMINT NOT NULL AUTO_INCREMENT,
            reactome_identifier VARCHAR(255) NOT NULL UNIQUE,
            reaction_name VARCHAR(255) NOT NULL,
            organism VARCHAR(255) NOT NULL,
            reaction_url VARCHAR(255),
            PRIMARY KEY(reaction_id)
        )
    '''

    return mysql_create_table(create_reaction_table_sql)

def create_uniprot_to_reaction_table():
   
    create_uniprot_to_reaction_table_sql = '''
        CREATE TABLE uniprot_to_reaction (
            uniprot_id MEDIUMINT NOT NULL,
            reaction_id MEDIUMINT NOT NULL,
            evidence VARCHAR(25),
            PRIMARY KEY(uniprot_id, reaction_id),
            FOREIGN KEY(uniprot_id) REFERENCES uniprot(uniprot_id),
            FOREIGN KEY(reaction_id) REFERENCES reaction(reaction_id)
        )
    '''

    return mysql_create_table(create_uniprot_to_reaction_table_sql)

def populate_uniprot_to_reactions():
    import pandas as pd

    uniprot_to_reaction = pd.read_csv(
        "/home/david/Desktop/uniprot_to_reaction.csv", 
        index_col=0)

    conn = connect_to_mysqldb()

    unique_uniprots = set(uniprot_to_reaction["acc"]) # relecant accs
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

    return mysql_insert_many(insert_uniprot_to_reaction_sql, rows, existing_conn=conn)

def create_compound_to_uniprot_table():

    create_compound_to_uniprot_table_sql = '''
        CREATE TABLE compound_to_uniprot (
            compound_id MEDIUMINT NOT NULL, 
            uniprot_id MEDIUMINT NOT NULL,
            confidence_score SMALLINT,
            PRIMARY KEY (uniprot_id, compound_id),
            FOREIGN KEY(compound_id) REFERENCES compounds(compound_id),
            FOREIGN KEY(uniprot_id) REFERENCES uniprot(uniprot_id),
            `above_950` bit(1) GENERATED ALWAYS AS (confidence_score>950) VIRTUAL NOT NULL,
            `above_900` bit(1) GENERATED ALWAYS AS (confidence_score>900) VIRTUAL NOT NULL,
            `above_850` bit(1) GENERATED ALWAYS AS (confidence_score>850) VIRTUAL NOT NULL,
            `above_800` bit(1) GENERATED ALWAYS AS (confidence_score>800) VIRTUAL NOT NULL,
            `above_750` bit(1) GENERATED ALWAYS AS (confidence_score>750) VIRTUAL NOT NULL,
            `above_700` bit(1) GENERATED ALWAYS AS (confidence_score>700) VIRTUAL NOT NULL,
            KEY `above_950_uniprot_id` (`above_950`, `uniprot_id`),
            KEY `above_900_uniprot_id` (`above_900`, `uniprot_id`),
            KEY `above_850_uniprot_id` (`above_850`, `uniprot_id`),
            KEY `above_800_uniprot_id` (`above_800`, `uniprot_id`),
            KEY `above_750_uniprot_id` (`above_750`, `uniprot_id`),
            KEY `above_700_uniprot_id` (`above_700`, `uniprot_id`)
        )

    '''

    return mysql_create_table(create_compound_to_uniprot_table_sql)

def populate_compound_to_uniprot_table(
    predictions_filename,
    min_confidence=650):
    import pandas as pd

    insert_compound_to_uniprot_sql = '''
        INSERT INTO compound_to_uniprot
            (compound_id, uniprot_id, confidence_score)
            VALUES (%s, %s, %s)
            ON DUPLICATE KEY UPDATE compound_id=compound_id
    '''

    print ("reading from", predictions_filename)
    predictions_df = pd.read_csv(predictions_filename, index_col=0)

    rows = ((compound_id, target_id, confidence_score)
        for compound_id, row in predictions_df.iterrows()
        for target_id, confidence_score in row.items()
        if confidence_score>=min_confidence)

    return mysql_insert_many(insert_compound_to_uniprot_sql, rows)

if __name__ == "__main__":
    # populate_targets_to_uniprot()

    # import pandas as pd
    # targets = pd.read_csv("models/target_ids.txt", 
    #     index_col=0, header=None, names=["uniprot"])
    # targets = tuple(targets["uniprot"])

    # add_uniprot_accs(targets)

    # create_compound_to_uniprot_table()
    for chunk_no in (1, ):
        predictions_filename = f"coconut_uniprot_predictions_chunk_{chunk_no}.csv.gz"
        populate_compound_to_uniprot_table(predictions_filename)