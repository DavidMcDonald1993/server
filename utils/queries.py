import sys
import os.path
sys.path.insert(1, 
    os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))

from collections import defaultdict

from utils.mysql_utils import mysql_query, sanitise_names, connect_to_mysqldb, mysql_create_table, mysql_execute

import urllib.parse as urlparse

import pandas as pd

def convert_to_dict(records, cols):
    return [
        {c:r for c, r in zip(cols, record)}
        for record in records
    ]

def get_all_targets_for_categories(
    categories=None, 
    existing_conn=None,
    as_dict=True):

    if categories is not None:
        if not isinstance(categories, str):
            if isinstance(categories, list) or isinstance(categories, set):
                categories = tuple(categories)
            assert isinstance(categories, tuple)
            if len(categories) == 1:
                categories = categories[0]

    query = f'''
       SELECT c.category_name AS 'target_category', 
       t.target_name AS `target_name`
       FROM categories AS c
       INNER JOIN category_members AS m ON (c.category_id=m.category_id) 
       INNER JOIN targets AS t ON (m.target_id=t.target_id)
       {f"WHERE c.category_name IN {categories}" 
        if isinstance(categories, tuple) 
        else f'WHERE c.category_name="{categories}"' 
            if isinstance(categories, str)
        else ""}
    '''

    records, cols = mysql_query(query, existing_conn=existing_conn, return_cols=True)
  
    if as_dict:
        records = convert_to_dict(records, cols)
  
    return records, cols

def get_uniprots_for_targets(targets, existing_conn=None, as_dict=False):
    if not isinstance(targets, str):
        if not isinstance(targets, tuple):
            targets = tuple(targets)
        assert isinstance(targets, tuple)
        print ("querying mySQL database for UNIPROTS associated with",
            len(targets), "targets")
        if len(targets) == 1:
            targets = targets[0]
    query = f'''
        SELECT DISTINCT t.target_name AS `target_name`, 
            u.acc AS `acc`, 
            u.protein AS `protein`, 
            u.gene AS `gene`, 
            tu.score AS `association_score`
        FROM targets AS t
        INNER JOIN targets_to_uniprot AS tu ON (t.target_id=tu.target_id)
        INNER JOIN uniprot AS u ON (tu.uniprot_id=u.uniprot_id)
        WHERE t.target_name {f"IN {targets}" if isinstance(targets, tuple)
            else f'="{targets}"'}
    '''
    return mysql_query(query, 
        existing_conn=existing_conn)

def get_targets_for_uniprot(accs, existing_conn=None, as_dict=True):
    if not isinstance(accs, str):
        if not isinstance(accs, tuple):
            accs = tuple(accs)
        assert isinstance(accs, tuple)
        if len(accs) == 1:
            accs = accs[0]

    query = f'''
        SELECT t.target_name AS `target_name`, 
            u.acc AS `acc`
        FROM targets AS t
        INNER JOIN targets_to_uniprot AS tu 
            ON (t.target_id=tu.target_id)
        INNER JOIN uniprot AS u 
            ON (tu.uniprot_id=u.uniprot_id)
        WHERE u.acc {f"IN {accs}" if isinstance(accs, tuple)
            else f'="{accs}"'}
    '''
    records, cols = mysql_query(query, 
        existing_conn=existing_conn, return_cols=True)
    if as_dict:
        records = convert_to_dict(records, cols)
    return records, cols


# def get_predicted_uniprots_for_compounds(
#     coconut_id, 
#     threshold=900,
#     filter_pa_pi=True,
#     as_dict=True,
#     existing_conn=None):

#     if isinstance(coconut_id, list) or isinstance(coconut_id, set):
#         coconut_id = tuple(coconut_id)

#     get_predicted_uniprots_sql = f'''
#     SELECT c.coconut_id AS `id`, 
#         u.acc AS `acc,
#         u.protein AS `protein`,
#         u.gene AS `gene`,
#         u.organism AS `organism`,
#         cu.confidence_score AS `confidence_score`
#     FROM compounds AS c
#     INNER JOIN compound_to_uniprot AS cu 
#         ON (c.compound_id=cu.compound_id)
#     INNER JOIN uniprot AS u 
#         ON (cu.uniprot_id=u.uniprot_id)
#     WHERE {f'c.coconut_id="{coconut_id}"' if isinstance(coconut_id, str)
#         else f"c.coconut_id IN {coconut_id}"}
#     AND cu.above_{threshold} = 1
#     ORDER BY confidence_score DESC
#     '''
#     return mysql_query(get_predicted_uniprots_sql, existing_conn=existing_conn)

# def get_inferred_uniprots_for_compounds(
#     coconut_id, 
#     threshold=750,
#     filter_pa_pi=True,
#     as_dict=True,
#     existing_conn=None):
#     assert filter_pa_pi

#     if isinstance(coconut_id, list) or isinstance(coconut_id, set):
#         coconut_id = tuple(coconut_id)

#     get_inferred_uniprots_sql = f'''
#     SELECT c.coconut_id AS `id`, 
#         u.acc AS `acc`,
#         u.protein AS `protein`,
#         u.gene AS `gene`,
#         u.organism AS `organism`,
#         FLOOR(MAX(a.confidence_score)) AS `max_confidence`,
#         GROUP_CONCAT(CONCAT("(", t.target_name, ",", a.confidence_score, ")") ORDER BY a.confidence_score DESC SEPARATOR '-') AS `targets`
#     FROM compounds AS c
#     INNER JOIN activities AS a ON (c.compound_id=a.compound_id)
#     INNER JOIN targets AS t ON (a.target_id=t.target_id)
#     INNER JOIN targets_to_uniprot AS tu ON (a.target_id=tu.target_id)
#     INNER JOIN uniprot AS u ON (tu.uniprot_id=u.uniprot_id)
#     WHERE {f"c.coconut_id='{coconut_id}'" if isinstance(coconut_id, str)
#         else f"c.coconut_id IN {coconut_id}"}
#     AND a.above_{threshold}=(1)
#     GROUP BY u.acc, u.protein, u.gene, u.organism
#     ORDER BY max_confidence DESC
#     '''
#     return mysql_query(get_inferred_uniprots_sql, existing_conn=existing_conn)

# def get_predicted_compounds_for_uniprots(
#     accs, 
#     threshold=900,
#     existing_conn=None):

#     if isinstance(accs, list) or isinstance(accs, set):
#         accs = tuple(accs)

#     get_predicted_compounds_sql = f'''
#     SELECT c.coconut_id AS `id`, 
#         c.image AS `image`,
#         c.name AS `name`, 
#         c.formula AS `formula`, 
#         c.smiles AS `smiles`,
#         u.acc AS `acc`,
#         cu.confidence_score AS `confidence`
#     FROM compounds AS c
#     INNER JOIN compound_to_uniprot AS cu 
#         ON (c.compound_id=cu.compound_id)
#     INNER JOIN uniprot AS u 
#         ON (cu.uniprot_id=u.uniprot_id)
#     WHERE {f'u.acc="{accs}"' if isinstance(accs, str)
#         else f"u.acc IN {accs}"}
#     {f"AND cu.above_{threshold}=(1)" if threshold > 900 else ""}
#     ORDER BY `confidence` DESC
#     '''
#     return mysql_query(get_predicted_compounds_sql, 
#         # return_cols=True,
#         existing_conn=existing_conn)

# def get_inferred_compounds_for_uniprots(
    # accs, 
    # threshold=750,
    # filter_pa_pi=True,
    # existing_conn=None):
    # assert filter_pa_pi

    # if isinstance(accs, list) or isinstance(accs, set):
    #     accs = tuple(accs)

    # get_inferred_compounds_sql = f'''
    # SELECT c.coconut_id AS `id`, 
    #     c.image AS `image`,
    #     c.name AS `name`, 
    #     c.formula AS `formula`, 
    #     c.smiles AS `smiles`,
    #     u.acc AS `acc`,
    #     FLOOR(MAX(a.confidence_score)) AS `max_confidence`,
    #     GROUP_CONCAT( CONCAT("(", t.target_name, ",", a.confidence_score, ")") ORDER BY a.confidence_score DESC SEPARATOR '-') AS `targets`
    # FROM compounds AS c
    # INNER JOIN activities AS a ON (c.compound_id=a.compound_id)
    # INNER JOIN targets AS t ON (a.target_id=t.target_id)
    # INNER JOIN targets_to_uniprot AS tu ON (a.target_id=tu.target_id)
    # INNER JOIN uniprot AS u ON (tu.uniprot_id=u.uniprot_id)
    # WHERE {f"u.acc='{accs}'" if isinstance(accs, str)
    #     else f"u.acc IN {accs}"}
    # AND a.above_{threshold}=(1)
    # GROUP BY c.compound_id
    # ORDER BY `max_confidence` DESC
    # '''
    # return mysql_query(get_inferred_compounds_sql, 
        # return_cols=True,
        # existing_conn=existing_conn)

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
        INNER JOIN uniprot_to_pathway AS up
            ON (p.pathway_id=up.pathway_id)
        INNER JOIN (
            SELECT DISTINCT uniprot_id
            FROM compound_to_uniprot AS cu
            UNION
            SELECT DISTINCT uniprot_id
            FROM targets_to_uniprot AS tu
        ) AS uniprots ON (up.uniprot_id=uniprots.uniprot_id)
    '''

    query = f'''
        SELECT DISTINCT p.pathway_name AS `pathway`, 
            p.organism
        FROM pathway AS p
        {active_filter if filter_actives else ""}
        {f"WHERE p.organism IN {organisms}" 
            if isinstance(organisms, tuple) else f'WHERE p.organism="{organisms}"' 
                if isinstance(organisms, str) else ""}
    '''
    return mysql_query(query, existing_conn=existing_conn)

def get_all_pathway_organisms(existing_conn=None):
    query = '''
        SELECT organism, COUNT(pathway_name) AS `num_pathways`
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
        INNER JOIN uniprot_to_reaction AS ur
            ON (r.reaction_id=ur.reaction_id)
        INNER JOIN (
            SELECT DISTINCT uniprot_id
            FROM compound_to_uniprot AS cu
            UNION
            SELECT DISTINCT uniprot_id
            FROM targets_to_uniprot AS tu
        ) AS uniprots ON (ur.uniprot_id=uniprots.uniprot_id)
    '''
    query = f'''
        SELECT DISTINCT r.reaction_name AS `reaction`, r.organism AS `organism`
        FROM reaction AS r
        {active_filter if filter_actives else ""}
        {f"WHERE r.organism IN {organisms}" 
            if isinstance(organisms, tuple) else f'WHERE r.organism="{organisms}"' 
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

def get_all_pathways_for_uniprots(
    accs,
    organism=None,
    existing_conn=None,
    as_dict=False,
    limit=None,
    ):
    if isinstance(accs, list) or isinstance(accs, set):
        accs = tuple(accs)
    else:
        assert isinstance(accs, str)

    query = f'''
    SELECT GROUP_CONCAT(DISTINCT(u.acc) SEPARATOR '-') AS `accs`, 
        COUNT(DISTINCT(u.acc)) AS `num_accs`,
        p.num_uniprot AS `total_accs`,
        COUNT(DISTINCT(u.acc)) / p.num_uniprot AS `acc_coverage`,      
        p.pathway_name AS `pathway`, p.organism AS `organism`, 
        p.pathway_url AS `url`
    FROM uniprot AS u
    INNER JOIN uniprot_to_pathway AS up ON (u.uniprot_id=up.uniprot_id)
    INNER JOIN pathway AS p ON (up.pathway_id=p.pathway_id)
    WHERE {f"u.acc IN {accs}" if isinstance(accs, tuple)
        else f'u.acc="{accs}"'}
    {f'AND p.organism="{organism}"' if organism is not None else ""}
    GROUP BY p.pathway_name, p.organism
    {f"LIMIT {limit}" if limit is not None else ""}
    '''
    records, cols = mysql_query(query, existing_conn=existing_conn, return_cols=True)

    # records = [
    #     (uniprots, uniprot_count, total_uniprots,
    #         coverage, pathway_name, urlparse.quote(pathway_name),
    #         organism, urlparse.quote(organism), url)
    #     for uniprots, uniprot_count, total_uniprots,
    #         coverage, pathway_name,
    #         organism, url in records
    # ]
    if as_dict:
        records = convert_to_dict(records, cols)

    return records, cols

def get_all_reactions_for_uniprots(
    accs,
    organism=None,
    existing_conn=None,
    as_dict=False,
    limit=None,
    ):
    if isinstance(accs, list) or isinstance(accs, set):
        accs = tuple(accs)
    else:
        assert isinstance(accs, str)
    query = f'''
    SELECT GROUP_CONCAT(DISTINCT(u.acc) SEPARATOR '-') AS `accs`, 
        COUNT(DISTINCT(u.acc)) AS `num_accs`,
        r.num_uniprot AS `total_accs`,
        COUNT(DISTINCT(u.acc)) / r.num_uniprot AS `acc_coverage`,   
        r.reaction_name AS `reacion`, r.organism AS `organism`, 
        r.reaction_url AS `url`
    FROM uniprot AS u
    INNER JOIN uniprot_to_reaction AS ur ON (u.uniprot_id=ur.uniprot_id)
    INNER JOIN reaction AS r ON (ur.reaction_id=r.reaction_id)
    WHERE {f"u.acc IN {accs}" if isinstance(accs, tuple)
        else f'u.acc="{accs}"'}
    {f'AND r.organism="{organism}"' if organism is not None else ""}
    GROUP BY r.reaction_name, r.organism
    {f"LIMIT {limit}" if limit is not None else ""}
    '''
    records, cols = mysql_query(query, existing_conn=existing_conn, return_cols=True)

    # records = [
    #     (uniprots, uniprot_count, total_uniprots,
    #         coverage, reaction_name, urlparse.quote(reaction_name,),
    #         organism, urlparse.quote(organism), url)
    #     for uniprots, uniprot_count, total_uniprots,
    #         coverage, reaction_name,
    #         organism, url in records
    # ]
    if as_dict:
        records = convert_to_dict(records, cols)

    return records, cols

def get_target_hits(
    targets, 
    threshold=750,
    min_targets_hit=None,
    as_dict=True,
    filter_pa_pi=True,
    limit=None
    ):
    assert filter_pa_pi
    assert isinstance(targets, list) or isinstance(targets, set)
    num_targets = len(targets)
    targets = tuple(targets)
    if len(targets)==1:
        targets = targets[0]
    assert isinstance(threshold, int)
    if min_targets_hit is not None:
        assert isinstance(min_targets_hit, int)

    query = f'''
        SELECT 
            c.image AS `image`, 
            c.coconut_id AS `id`, 
            c.name AS `name`, 
            c.formula AS `formula`, 
            c.smiles AS `smiles`,
            t.target_name AS `target_name`,
            a.confidence_score AS `confidence_score`
        FROM compounds AS `c`
        INNER JOIN activities AS `a`
            ON (c.compound_id=a.compound_id)
        INNER JOIN targets AS `t`
            ON (a.target_id=t.target_id)
        WHERE {f"t.target_name IN {targets}" if isinstance(targets, tuple)
            else f't.target_name="{targets}"'}
        AND a.above_{threshold}=(1)
        ORDER BY confidence_score DESC
    '''

    if isinstance(targets, tuple):
        query = f'''
        SELECT id, image, name, formula, smiles,
        COUNT(target_name) AS `num_targets_hit`,
        FLOOR(AVG(confidence_score)) AS `mean_confidence`,
        JSON_ARRAYAGG(
        JSON_OBJECT(
            'target_name', target_name,
            'confidence_score', confidence_score
        )
        ) AS all_confidences
        FROM ({query}) AS predictions
        GROUP BY id, image, name, formula, smiles
        {f"HAVING num_targets_hit>={min_targets_hit}" if min_targets_hit is not None else ""}
        ORDER BY num_targets_hit DESC, mean_confidence DESC
        '''

    records, cols = mysql_query(query, return_cols=True)

    if as_dict:
        records = convert_to_dict(records, cols)

    return records, cols

def get_pathway_hits(
    pathways,
    organism, 
    threshold=750, 
    filter_pa_pi=True, 
    min_target_coverage=None,
    min_pathways_hit=None,
    as_dict=True,
    limit=None,
    existing_conn=None):
    assert filter_pa_pi

    if isinstance(pathways, list) or isinstance(pathways, set):
        pathways = tuple(pathways)
        if len(pathways) == 1:
            pathways = pathways[0]
    else:
        assert isinstance(pathways, str)

    if existing_conn is None:
        existing_conn = connect_to_mysqldb()

    get_pathway_id_sql = f'''
        SELECT pathway_id
        FROM pathway AS p
        WHERE {f"p.pathway_name IN {pathways}" if isinstance(pathways, tuple)
            else f'p.pathway_name="{pathways}"'}
        AND p.organism="{organism}"
    '''

    records = mysql_query(get_pathway_id_sql, existing_conn=existing_conn)

    if len(records) == 1:
        ids = records[0][0]
    else:
        ids = tuple([record[0] for record in records])

    predicted_query = f'''
        SELECT cu.compound_id AS `predicted_compound_id`, 
            cu.uniprot_id AS `predicted_uniprot_id`, 
            cu.confidence_score AS `predicted_confidence_score`, 
            up.pathway_id AS `predicted_pathway_id`,
            up.evidence AS `predicted_pathway_evidence`
        FROM compound_to_uniprot AS cu 
        INNER JOIN uniprot_to_pathway AS up
            ON (cu.uniprot_id=up.uniprot_id)
        WHERE {f'up.pathway_id={ids}'
            if isinstance(ids, int) else
            f"up.pathway_id IN {ids}"}
        AND cu.above_{threshold}=(1)
    '''

    inferred_query = f'''
        SELECT a.compound_id AS `inferred_compound_id`, 
            tu.uniprot_id AS `inferred_uniprot_id`, 
            FLOOR(MAX(a.confidence_score)) AS `inferred_confidence_score`, 
            GROUP_CONCAT(t.target_name) AS `inferred_relevant_targets`,
            up.pathway_id AS `inferred_pathway_id`,
            up.evidence AS `inferred_pathway_evidence`
        FROM activities AS a
        INNER JOIN targets AS t
            ON (a.target_id=t.target_id)
        INNER JOIN targets_to_uniprot as tu 
            on (a.target_id=tu.target_id) 
        INNER JOIN uniprot_to_pathway AS up
            ON (tu.uniprot_id=up.uniprot_id)
        WHERE {f'up.pathway_id={ids}'
            if isinstance(ids, int) else
            f"up.pathway_id IN {ids}"}
        AND a.above_{threshold}=(1) 
        GROUP BY inferred_compound_id, inferred_uniprot_id, inferred_pathway_id, inferred_pathway_evidence
    '''

    create_combined_sql = f'''
        SELECT
            CASE
                WHEN predictions.predicted_compound_id IS NOT NULL 
                    THEN  predictions.predicted_compound_id
                ELSE predictions.inferred_compound_id
            END AS `compound_id`,
            CASE 
                WHEN predictions.predicted_pathway_id  IS NOT NULL
                    THEN predictions.predicted_pathway_id 
                ELSE predictions.inferred_pathway_id 
            END AS `pathway_id`,
            CASE 
                WHEN predictions.predicted_uniprot_id IS NOT NULL THEN predictions.predicted_uniprot_id
                ELSE predictions.inferred_uniprot_id
            END AS `uniprot_id`,
            CASE 
                WHEN predictions.predicted_pathway_evidence IS NOT NULL THEN predictions.predicted_pathway_evidence
                ELSE predictions.inferred_pathway_evidence
            END AS `evidence`,
            CASE 
                WHEN predictions.predicted_confidence_score IS NOT NULL AND 
                    predictions.inferred_confidence_score IS NOT NULL THEN 1000
                WHEN predictions.predicted_confidence_score IS NOT NULL
                    THEN predictions.predicted_confidence_score
                ELSE predictions.inferred_confidence_score
            END AS `confidence_score`,
            CASE 
                WHEN predictions.predicted_confidence_score IS NOT NULL AND 
                    predictions.inferred_confidence_score IS NOT NULL THEN "CONFIRMED"
                WHEN predictions.predicted_confidence_score IS NOT NULL
                    THEN "PREDICTED"
                ELSE "INFERRED"
            END AS `confidence_type`,
            predictions.inferred_relevant_targets
        FROM (
            WITH predicted AS ({predicted_query}), inferred AS ({inferred_query})
            SELECT *
            FROM predicted
            LEFT JOIN 
            inferred
                ON (predicted.predicted_uniprot_id=inferred.inferred_uniprot_id
                    AND predicted.predicted_compound_id=inferred.inferred_compound_id
                    AND predicted.predicted_pathway_id=inferred.inferred_pathway_id)
            UNION ALL
            SELECT * 
            FROM predicted
            RIGHT JOIN
            inferred
                ON (predicted.predicted_uniprot_id=inferred.inferred_uniprot_id
                    AND predicted.predicted_compound_id=inferred.inferred_compound_id
                    AND predicted.predicted_pathway_id=inferred.inferred_pathway_id)
            WHERE predicted.predicted_uniprot_id IS NULL
        ) AS predictions
    '''

    # combined_table_name = "temp_combined_pathway"
    # create_temp_combined_table_sql = f'''
    # CREATE TEMPORARY TABLE {combined_table_name} (
    #     compound_id MEDIUMINT, 
    #     pathway_id MEDIUMINT,
    #     uniprot_id MEDIUMINT,
    #     evidence VARCHAR(5),
    #     confidence_score SMALLINT,
    #     confidence_type VARCHAR(10),
    #     PRIMARY KEY(compound_id, uniprot_id, pathway_id),
    #     KEY (compound_id),
    #     KEY (uniprot_id),
    #     KEY (pathway_id)
    # )
    # {create_combined_sql}
    # '''
    # mysql_create_table(create_temp_combined_table_sql, existing_conn=existing_conn)

    query = f'''
        SELECT 
            c.image AS `image`, 
            c.coconut_id AS `id`, 
            c.name AS `name`, 
            c.formula AS `formula`, 
            c.smiles AS `smiles`,
            p.pathway_name AS `pathway_name`,
            p.organism AS `pathway_organism`,
        COUNT(u.acc) AS num_targets_hit,
        p.num_uniprot AS total_targets,
        CAST(COUNT(u.acc) / p.num_uniprot AS CHAR) AS `coverage`,
        JSON_ARRAYAGG(
            JSON_OBJECT(
                'target_acc', u.acc,
                'target_protein', u.protein,
                'target_gene', u.gene,
                'prediction_confidence', predictions.confidence_score,
                'annotation_type', predictions.confidence_type,
                'evidence', predictions.evidence,
                'relevant_targets', predictions.inferred_relevant_targets
            )
        ) AS associated_predicted_targets,
        p.pathway_url AS url
        FROM ({create_combined_sql}) AS predictions
        INNER JOIN compounds AS c
            ON (c.compound_id=predictions.compound_id)
        INNER JOIN uniprot AS u
            ON (predictions.uniprot_id=u.uniprot_id)
        INNER JOIN pathway AS p
            ON (predictions.pathway_id=p.pathway_id)
        GROUP BY id, image, name, formula, smiles, pathway_name
        {f"HAVING coverage>={min_target_coverage}" if min_target_coverage is not None else ""}
        ORDER BY coverage DESC
    '''

    # records = mysql_query(query, existing_conn=existing_conn)

    # for record in records[:5]:
    #     print (record)

    # mysql_execute("DROP TABLE predicted", existing_conn=existing_conn)
    # mysql_execute("DROP TABLE inferred", existing_conn=existing_conn)
    # raise Exception

    #  GROUP_CONCAT( 
    #         CONCAT( "{{pathway_name:", predictions.pathway_name, ", target_coverage:", predictions.coverage, ", accs:", predictions.accs, "}}" )
    #         ORDER BY predictions.coverage DESC
    #         SEPARATOR ","
        # ) AS `summary`

    if isinstance(pathways, tuple):
        query = f'''
            SELECT 
                predictions.`image`, 
                predictions.`id`, 
                predictions.`name`, 
                predictions.`formula`, 
                predictions.`smiles`,
            GROUP_CONCAT(predictions.pathway_name) AS `pathways_hit`,
            COUNT(predictions.pathway_name) AS `num_pathways_hit`,
            CAST(ROUND(AVG(predictions.coverage), 3) AS CHAR) AS `expected_coverage`,
            JSON_ARRAYAGG(
                JSON_OBJECT(
                    'pathway_name', predictions.pathway_name,
                    'pathway_organism', predictions.pathway_organism,
                    'num_targets_hit', predictions.num_targets_hit,
                    'total_targets', predictions.total_targets,
                    'target_coverage', predictions.coverage,
                    'associated_predicted_targets', predictions.associated_predicted_targets,
                    'url', predictions.url
                )
            ) AS summary
            FROM ({query}) AS `predictions`
            GROUP BY id, image, name, formula, smiles
            {f"HAVING num_pathways_hit>={min_pathways_hit}" if min_pathways_hit is not None else ""}
            ORDER BY num_pathways_hit DESC, expected_coverage ASC
        '''

    records, cols = mysql_query(query, return_cols=True, existing_conn=existing_conn)

    if as_dict:
        records = [
            {c:r for c, r in zip(cols, record)}
            for record in records
        ]

    return records, cols

def get_reaction_hits(
    reactions,
    organism,
    threshold=750,
    filter_pa_pi=True, 
    min_target_coverage=None,
    min_reactions_hit=None,
    as_dict=True,
    limit=None,
    existing_conn=None):
    assert filter_pa_pi

    if isinstance(reactions, list) or isinstance(reactions, set):
        reactions = tuple(reactions)
        if len(reactions) == 1:
            reactions = reactions[0]
    else:
        assert isinstance(reactions, str)

    if existing_conn is None:
        existing_conn = connect_to_mysqldb()

    get_reaction_id_sql = f'''
        SELECT reaction_id
        FROM reaction AS r
        WHERE {f"r.reaction_name IN {reactions}" if isinstance(reactions, tuple)
            else f'r.reaction_name="{reactions}"'}
        AND r.organism="{organism}"
    '''

    records = mysql_query(get_reaction_id_sql, existing_conn=existing_conn)

    if len(records) == 1:
        ids = records[0][0]
    else:
        ids = tuple([record[0] for record in records])

    predicted_query = f'''
        SELECT cu.compound_id AS `predicted_compound_id`, 
            cu.uniprot_id AS `predicted_uniprot_id`, 
            cu.confidence_score AS `predicted_confidence_score`, 
            ur.reaction_id AS `predicted_reaction_id`,
            ur.evidence AS `predicted_reaction_evidence`
        FROM compound_to_uniprot AS cu 
        INNER JOIN uniprot_to_reaction AS ur
            ON (cu.uniprot_id=ur.uniprot_id)
        WHERE {f'ur.reaction_id={ids}'
            if isinstance(ids, int) else
            f"ur.reaction_id IN {ids}"}
        AND cu.above_{threshold}=(1)
    '''

    inferred_query = f'''
        SELECT a.compound_id AS `inferred_compound_id`, 
            tu.uniprot_id AS `inferred_uniprot_id`, 
            FLOOR(MAX(a.confidence_score)) AS `inferred_confidence_score`, 
            GROUP_CONCAT(t.target_name) AS `inferred_relevant_targets`,
            ur.reaction_id AS `inferred_reaction_id`,
            ur.evidence AS `inferred_reaction_evidence`
        FROM activities AS a
        INNER JOIN targets_to_uniprot as tu 
            on (a.target_id=tu.target_id) 
        INNER JOIN targets AS t
            ON (a.target_id=t.target_id)
        INNER JOIN uniprot_to_reaction AS ur
            ON (tu.uniprot_id=ur.uniprot_id)
        WHERE {f'ur.reaction_id={ids}'
            if isinstance(ids, int) else
            f"ur.reaction_id IN {ids}"}
        AND a.above_{threshold}=(1) 
        GROUP BY inferred_compound_id, inferred_uniprot_id, inferred_reaction_id, inferred_reaction_evidence
    '''

    create_combined_sql = f'''
        SELECT 
            CASE
                WHEN predictions.predicted_compound_id IS NOT NULL 
                    THEN  predictions.predicted_compound_id
                ELSE predictions.inferred_compound_id
            END AS `compound_id`,
            CASE 
                WHEN predictions.predicted_reaction_id IS NOT NULL
                    THEN predictions.predicted_reaction_id 
                ELSE predictions.inferred_reaction_id 
            END AS `reaction_id`,
            CASE 
                WHEN predictions.predicted_uniprot_id IS NOT NULL THEN predictions.predicted_uniprot_id
                ELSE predictions.inferred_uniprot_id
            END AS `uniprot_id`,
            CASE 
                WHEN predictions.predicted_reaction_evidence IS NOT NULL THEN predictions.predicted_reaction_evidence
                ELSE predictions.inferred_reaction_evidence
            END AS `evidence`,
            CASE 
                WHEN predictions.predicted_confidence_score IS NOT NULL AND 
                    predictions.inferred_confidence_score IS NOT NULL THEN 1000
                WHEN predictions.predicted_confidence_score IS NOT NULL
                    THEN predictions.predicted_confidence_score
                ELSE predictions.inferred_confidence_score
            END AS `confidence_score`,
            CASE 
                WHEN predictions.predicted_confidence_score IS NOT NULL AND 
                    predictions.inferred_confidence_score IS NOT NULL THEN "CONFIRMED"
                WHEN predictions.predicted_confidence_score IS NOT NULL
                    THEN "PREDICTED"
                ELSE "INFERRED"
            END AS `confidence_type`,
            predictions.inferred_relevant_targets
        FROM (
            WITH predicted AS ({predicted_query}), inferred AS ({inferred_query})
            SELECT *
            FROM predicted
            LEFT JOIN 
            inferred
                ON (predicted.predicted_uniprot_id=inferred.inferred_uniprot_id
                    AND predicted.predicted_compound_id=inferred.inferred_compound_id
                    AND predicted.predicted_reaction_id=inferred.inferred_reaction_id)
            UNION ALL
            SELECT * 
            FROM  predicted
            RIGHT JOIN 
            inferred
                ON (predicted.predicted_uniprot_id=inferred.inferred_uniprot_id
                    AND predicted.predicted_compound_id=inferred.inferred_compound_id
                    AND predicted.predicted_reaction_id=inferred.inferred_reaction_id)
            WHERE predicted.predicted_uniprot_id IS NULL
        ) AS predictions
    '''
    # combined_table_name = "temp_combined_reaction"
    # create_temp_combined_table_sql = f'''
    # CREATE TEMPORARY TABLE {combined_table_name} (
    #     compound_id MEDIUMINT, 
    #     reaction_id MEDIUMINT,
    #     uniprot_id MEDIUMINT,
    #     evidence VARCHAR(5),
    #     confidence_score SMALLINT,
    #     confidence_type VARCHAR(10),
    #     KEY (compound_id),
    #     KEY (uniprot_id),
    #     KEY (reaction_id)
    # )
    # {create_combined_sql}
    # '''
    # mysql_create_table(create_temp_combined_table_sql, existing_conn=existing_conn)

    query = f'''
        SELECT 
            c.image AS `image`, 
            c.coconut_id AS `id`, 
            c.name AS `name`, 
            c.formula AS `formula`, 
            c.smiles AS `smiles`,
            r.reaction_name,
            r.organism AS `reaction_organism`,
            COUNT(u.acc) AS num_targets_hit,
            r.num_uniprot AS total_targets,
            CAST(COUNT(u.acc) / r.num_uniprot AS CHAR) AS `coverage`,
            JSON_ARRAYAGG(
                JSON_OBJECT(
                    'target_acc', u.acc,
                    'target_protein', u.protein,
                    'target_gene', u.gene,
                    'prediction_confidence', predictions.confidence_score,
                    'annotation_type', predictions.confidence_type,
                    'evidence', predictions.evidence,
                    'relevant_targets', predictions.inferred_relevant_targets
                )
            ) AS associated_predicted_targets,
            r.reaction_url AS url
        FROM 
        ({create_combined_sql}) AS predictions
        INNER JOIN compounds AS c
            ON (c.compound_id=predictions.compound_id)
        INNER JOIN uniprot AS u
            ON (predictions.uniprot_id=u.uniprot_id)
        INNER JOIN reaction AS r
            ON (predictions.reaction_id=r.reaction_id)
        GROUP BY id, reaction_name
        {f"HAVING coverage>={min_target_coverage}" if min_target_coverage is not None else ""}
        ORDER BY coverage DESC
    '''

    if isinstance(reactions, tuple):
        query = f'''
            SELECT 
                predictions.image AS `image`, 
                predictions.id AS `id`, 
                predictions.name AS `name`, 
                predictions.formula AS `formula`, 
                predictions.smiles AS `smiles`,
                GROUP_CONCAT(predictions.reaction_name) AS `reactions_hit`,
                COUNT(predictions.reaction_name) AS `num_reactions_hit`,
                CAST(ROUND(AVG(predictions.coverage), 3) AS CHAR) AS `expected_coverage`,
                JSON_ARRAYAGG(
                    JSON_OBJECT(
                        'reaction_name', predictions.reaction_name,
                        'reaction_organism', predictions.reaction_organism,
                        'target_coverage', predictions.coverage,
                        'num_targets_hit', predictions.num_targets_hit,
                        'total_targets', predictions.total_targets,
                        'associated_predicted_targets', predictions.associated_predicted_targets,
                        'url', predictions.url
                    )
                ) AS summary
            FROM ({query}) AS `predictions`
            GROUP BY id, image, name, formula, smiles
            {f"HAVING num_reactions_hit>={min_reactions_hit}" if min_reactions_hit is not None else ""}
            ORDER BY num_reactions_hit DESC, expected_coverage ASC
        '''

    records, cols = mysql_query(query, return_cols=True, existing_conn=existing_conn)

    if as_dict:
        records = [
            {c:r for c, r in zip(cols, record)}
            for record in records
        ]

    return records, cols

def get_all_kingdoms():

    query = '''
    SELECT kingdom_name
    FROM kingdom
    '''

    return [
        record[0] for record in mysql_query(query)
    ]

def get_all_species(species_group=None, limit=None):

    if species_group is not None:
        if isinstance(species_group, list) or isinstance(species_group, set):
            species_group = tuple(species_group)
            if len(species_group) == 1:
                species_group = species_group[0]

    query = f'''
    SELECT species_name
    FROM species
    {f"WHERE species_group IN {species_group}" if isinstance(species_group, tuple)
        else f'WHERE species_group="{species_group}"' if isinstance(species_group, str)
        else ""}
    {f"LIMIT {limit}" if limit is not None else ""}
    '''

    return [
        record[0] for record in mysql_query(query)
    ]

def get_all_species_groups(limit=None):

    query = f'''
    SELECT DISTINCT species_group
    FROM species
    {f"LIMIT {limit}" if limit is not None else ""}
    '''

    return [
        record[0] for record in mysql_query(query)
    ]

def get_info_for_multiple_compounds(
    compound_ids=None, 
    name_like=None,
    formula_like=None,
    smiles_like=None,
    kingdom_name=None,
    species_group=None,
    species_name=None,
    columns=("coconut_id", "name", "formula", "smiles", "kingdom_name", "species_name", ),
    as_dict=False,
    limit=None):

    if compound_ids is not None:
        if not isinstance(compound_ids, str):
            if isinstance(compound_ids, list) or isinstance(compound_ids, set):
                compound_ids = tuple(compound_ids)
            assert isinstance(compound_ids, tuple)
            if len(compound_ids) == 1:
                compound_ids = compound_ids[0]

    if isinstance(kingdom_name, list) or isinstance(kingdom_name, set):
        kingdom_name = tuple(kingdom_name)
        if len(kingdom_name) == 1:
            kingdom_name = kingdom_name[0]

    if isinstance(species_group, list) or isinstance(species_group, set):
        species_group = tuple(species_group)
        if len(species_group) == 1:
            species_group = species_group[0]

    if isinstance(species_name, list) or isinstance(species_name, set):
        species_name = tuple(species_name)
        if len(species_name) == 1:
            species_name = species_name[0]

    conditions = [
        f"coconut_id IN {compound_ids}" if isinstance(compound_ids, tuple)
            else f'coconut_id="{compound_ids}"' if isinstance(compound_ids, str) else None,
        f'name LIKE "%{name_like}%"' if name_like is not None else None,
        f'formula LIKE "%{formula_like}%"' if formula_like is not None else None,
        f'smiles LIKE "%{smiles_like}%"' if smiles_like is not None else None,
        f'kingdom_name = "{kingdom_name}"' if isinstance(kingdom_name, str)
            else f"kingdom_name IN {kingdom_name}" if isinstance(kingdom_name, tuple) else None,
        f'species_group = "{species_group}"' if isinstance(species_group, str)
            else f"species_group IN {species_group}" if isinstance(species_group, tuple) else None,
        f'species_name = "{species_name}"' if isinstance(species_name, str)
            else f"species_name IN {species_name}" if isinstance(species_name, tuple) else None,
    ]
    conditions = " AND ".join(filter(lambda x: x, conditions))

    if "kingdom_name" or "species_name" in columns:
        group_by = ",".join((c for c in columns
            if c not in {"kingdom_name", "species_name"}))
        
        columns = [
            f"GROUP_CONCAT(DISTINCT {c}) AS `all_{c}s`"
                if c in {"species_name", "kingdom_name"}
                else c 
            for c  in columns
        ]

    else: 
        group_by = None

    compound_query = f'''
        SELECT {(", ".join(columns))}
        FROM compounds AS c
        LEFT JOIN (
            SELECT compound_id AS cid, kingdom_name
            FROM kingdom AS k
            INNER JOIN compound_to_kingdom AS ck
                ON (k.kingdom_id=ck.kingdom_id)
        ) AS kingdom_link 
            ON (c.compound_id=kingdom_link.cid)
        LEFT JOIN (
            SELECT compound_id AS cid, species_name, species_group
            FROM species AS s
            INNER JOIN compound_to_species AS cs
                ON (cs.species_id=s.species_id)
        ) AS species_link
            ON (c.compound_id=species_link.cid)
        {f"WHERE {conditions}" if conditions != "" else ""}
        {f"GROUP BY {group_by}" if group_by is not None else ""}
        {f"LIMIT {limit}" if limit is not None else ""}
    '''

    records, cols = mysql_query(compound_query, return_cols=True)
    if as_dict:
        records = convert_to_dict(records, cols)
    return records, cols

def get_categories():
    query = '''
    SELECT category_name
    FROM categories
    ORDER BY category_name ASC
    '''
    return [record[0] for record in mysql_query(query)]

def get_all_activities_for_compound(
    coconut_id, 
    threshold=750,
    filter_pa_pi=True):
    assert filter_pa_pi
    print ("getting all activities for compound", coconut_id)

    all_targets_query = f'''
    SELECT cat.category_name, t.target_name, a.Pa, a.Pi, a.confidence_score
    FROM categories AS cat
    INNER JOIN category_members AS cm ON (cat.category_id=cm.category_id)
    INNER JOIN targets AS t ON (cm.target_id=t.target_id) 
    INNER JOIN activities AS a ON (t.target_id=a.target_id)
    INNER JOIN compounds AS c ON (a.compound_id=c.compound_id)
    WHERE c.coconut_id="{coconut_id}"
    {f"AND a.above_{threshold}=(1)" if threshold is not None and threshold > 0 and filter_pa_pi else ""}
    ORDER BY confidence_score DESC
    '''
    compound_hits = mysql_query(all_targets_query)

    compound_hits = [ # encode target name
        (category_name, target_name, urlparse.quote(target_name),
            pa, pi, confidence_score)
            for category_name, target_name, pa, pi, confidence_score
            in compound_hits
    ]
    
    category_activities = defaultdict(list)
    prev_hit = None # duplicates will appear next to each other
    for category, *hit in compound_hits:
        hit = tuple(hit)
        if prev_hit is None or prev_hit != hit:
            category_activities["ALL_TARGETS"].append(hit)
            prev_hit = hit
        category_activities[category].append(hit)
   
    return dict(category_activities)

def get_drugs_for_uniprots(accs, as_dict=True):
    if isinstance(accs, list) or isinstance(accs, set):
        accs = tuple(accs)
    else:
        assert isinstance(accs, str)
    uniprot_to_drug_sql = f'''
        SELECT ud.uniprot_id,
        ud.drug_id,
        ud.activity,
        ud.reference
        FROM uniprot AS u
        INNER JOIN uniprot_to_drug AS ud
            ON (u.uniprot_id=ud.uniprot_id)
        WHERE {f"u.acc IN {accs}" if isinstance(accs, tuple)
            else f'u.acc="{accs}"'}
    '''

    drug_ids = {record[1] for record in mysql_query(uniprot_to_drug_sql)}
    
    drug_ids = tuple(drug_ids)
    if len(drug_ids) == 1:
        drug_ids = drug_ids[0]

    if len(drug_ids) == 0:
        return None, None

    # disease_sql = f'''
    #     SELECT udis.uniprot_id, 
    #     drd.drug_id, 
    #     drd.clinical_status, 
    #     disease.disease_name
    #     FROM drug_to_disease as `drd`
    #     INNER JOIN uniprot_to_disease as `udis`
    #         ON (udis.disease_id=drd.disease_id)
    #     INNER JOIN disease
    #         ON (udis.disease_id=disease.disease_id)
    #     INNER JOIN uniprot AS u
    #         ON (u.uniprot_id=udis.uniprot_id)
    #     WHERE {f"u.acc IN {accs}" if isinstance(accs, tuple)
    #         else f'u.acc="{accs}"'}
    # '''

    # get all ossciated diseases for drugs
    disease_sql = f'''
        SELECT 
        u.acc,
        drd.drug_id, 
        drd.clinical_status, 
        disease.disease_name,
        disease.icd
        FROM drug_to_disease as `drd`
        INNER JOIN uniprot_to_disease as `udis`
            ON (udis.disease_id=drd.disease_id)
        INNER JOIN uniprot_to_drug AS ud
            ON (ud.uniprot_id=udis.uniprot_id
                AND ud.drug_id=drd.drug_id)
        INNER JOIN disease
            ON (udis.disease_id=disease.disease_id)
        INNER JOIN uniprot AS u
            ON (u.uniprot_id=udis.uniprot_id)
        WHERE {f"drd.drug_id IN {drug_ids}" if isinstance(drug_ids, tuple)
            else f'drd.drug_id="{drug_ids}"'}
        GROUP BY drug_id, disease_name
    '''

    # disease_sql = f'''
    #     SELECT CONCAT("(", GROUP_CONCAT(uniprot_id), ")") AS targets,
    #     drug_id, 
    #     clinical_status, 
    #     disease_name,
    #     icd
    #     FROM ({disease_sql}) AS t
    #     GROUP BY drug_id, disease_name
    # '''

    # records = mysql_query(disease_sql)
    # for record in records[:25]:
    #     print (record)
    # print (len(records))
    # raise Exception

    sql = f'''
    SELECT d.drug_name AS `drug_name`, 
        d.inchi AS `inchi`, 
        d.canonical_smiles AS `smiles`, 
        d.drug_type AS `type`, 
        d.drug_class AS `class`, 
        d.company AS `company`, 
        u.acc AS `acc`, 
        to_drug.activity AS `activity`, 
        to_drug.reference AS `reference`,
        CASE WHEN disease_link.disease_name IS NOT NULL
        THEN
            JSON_ARRAYAGG(
                JSON_OBJECT(
                    'disease_name', disease_link.disease_name,
                    'icd', disease_link.icd,
                    'clinical_status', disease_link.clinical_status,
                    'target', 
                    CASE
                        WHEN u.acc = disease_link.acc THEN CONCAT(u.acc, "(THIS TARGET)")
                        ELSE disease_link.acc
                    END
                )
            ) 
        ELSE NULL
        END AS all_diseases_for_drug
    FROM ({uniprot_to_drug_sql}) AS to_drug
    INNER JOIN uniprot AS u
        ON (to_drug.uniprot_id=u.uniprot_id)
    INNER JOIN drug AS d
        ON (to_drug.drug_id=d.drug_id)
    LEFT JOIN 
        ({disease_sql}) AS disease_link 
        ON (to_drug.drug_id=disease_link.drug_id)
    GROUP BY acc, drug_name
    '''
    records, cols = mysql_query(sql, return_cols=True)
    if as_dict:
        records = [
            {c:r for c, r in zip(cols, record)}
            for record in records
        ]
    return records, cols

def get_diseases_for_uniprots(accs, as_dict=True):
    if isinstance(accs, list) or isinstance(accs, set):
        accs = tuple(accs)
    else:
        assert isinstance(accs, str)
    sql = f'''
    SELECT 
        u.acc AS `acc`, 
        JSON_ARRAYAGG(
            JSON_OBJECT(
                "disease_name", d.disease_name,
                "icd", d.icd,
                "clinical_status", ud.clinical_status
            )
        ) AS `all_associated_diseases_for_target`
    FROM uniprot AS `u`
    INNER JOIN uniprot_to_disease as `ud`
        ON (u.uniprot_id=ud.uniprot_id)
    INNER JOIN disease as `d`
        ON (ud.disease_id=d.disease_id)
    WHERE {f"u.acc IN {accs}" if isinstance(accs, tuple)
        else f'u.acc="{accs}"'}
    GROUP BY acc
    '''
    records, cols = mysql_query(sql, return_cols=True)
    if as_dict:
        records = [
            {c:r for c, r in zip(cols, record)}
            for record in records
        ]
    return records, cols

def get_diseases_for_drugs(drugs, as_dict=True):
    if isinstance(drugs, list) or isinstance(drugs, set):
        drugs = tuple(drugs)
    else:
        assert isinstance(drugs, str)
    sql = f'''
    SELECT disease.disease_name AS `disease`, 
        disease.icd AS `icd`,
        d.drug_name AS `drug_name`, 
        drd.clinical_status AS `clinical_status`
    FROM drug AS `d`
    INNER JOIN drug_to_disease as `drd`
        ON (d.drug_id=drd.drug_id)
    INNER JOIN disease
        ON (drd.disease_id=disease.disease_id)
    WHERE {f"d.drug_name IN {drugs}" if isinstance(drugs, tuple)
        else f'd.drug_name="{drugs}"'}
    '''
    records, cols = mysql_query(sql, return_cols=True)
    if as_dict:
        records = [
            {c:r for c, r in zip(cols, record)}
            for record in records
        ]
    return records, cols

def get_protein_gene_from_acc(accs, as_dict=False):
    if isinstance(accs, list) or isinstance(accs, set):
        accs = tuple(accs)
    else:
        assert isinstance(accs, str)
    query = f'''
    SELECT acc, protein, gene
    FROM uniprot AS `u`
    WHERE {f"u.acc IN {accs}" if isinstance(accs, tuple)
        else f'u.acc="{accs}"'}
    '''
    records, cols = mysql_query(query, return_cols=True)
    if as_dict:
        records = convert_to_dict(records, cols)
    return records

def get_all_uniprot(as_dict=True, filter_valid=True):

    filter_valid_sql = '''
    WHERE u.uniprot_id IN
    (
        SELECT uniprot_id FROM compound_to_uniprot
        UNION 
        SELECT uniprot_id FROM targets_to_uniprot
        UNION
        SELECT uniprot_id FROM uniprot_to_drug
        UNION
        SELECT uniprot_id FROM uniprot_to_disease
    )
    '''

    query = f'''
    SELECT u.acc AS `acc`, 
        u.protein AS `protein`, 
        u.gene AS `gene`, 
        u.organism_common,
        u.organism_scientific,
        u.organism_synonym
    FROM uniprot as `u`
    {filter_valid_sql if filter_valid else ""}
    '''
    records, cols =  mysql_query(query, return_cols=True)
    if as_dict:
        records = convert_to_dict(records, cols)
    return records, cols

def get_combined_compound_confidences_for_uniprots(
    accs, 
    as_dict=True,
    existing_conn=None,
    threshold=750):

    if isinstance(accs, list) or isinstance(accs, set):
        accs = tuple(accs)
        if len(accs) == 1:
            accs = accs[0]
    else:
        assert isinstance(accs, str)

    if existing_conn is None:
        existing_conn = connect_to_mysqldb()

    predicted_query = f'''
        SELECT cu.compound_id AS `predicted_compound_id`, 
            u.acc AS `predicted_uniprot_acc`, 
            cu.confidence_score AS `predicted_confidence_score`
        FROM compound_to_uniprot AS cu 
        INNER JOIN uniprot AS u
            ON (cu.uniprot_id=u.uniprot_id)
        WHERE {f'u.acc="{accs}"'
            if isinstance(accs, str) else
            f"u.acc IN {accs}"}
        AND cu.above_{threshold}=(1)
    '''

    # records = mysql_query(predicted_query, existing_conn=existing_conn)

    # for record in records[:5]:
    #     print (record)
    # print (len(records))
    # existing_conn.close()
    # raise Exception

    inferred_query = f'''
        SELECT a.compound_id AS `inferred_compound_id`, 
            u.acc AS `inferred_uniprot_acc`, 
            MAX(a.confidence_score) AS `inferred_confidence_score`,
            JSON_ARRAYAGG(
                JSON_OBJECT(
                    'target_name', target_name,
                    'confidence_score', confidence_score
                )
            ) AS all_targets
        FROM activities AS a
        INNER JOIN targets_to_uniprot as tu 
            on (a.target_id=tu.target_id) 
        INNER JOIN targets AS t 
            ON (t.target_id=a.target_id)
        INNER JOIN uniprot AS u
            ON (tu.uniprot_id=u.uniprot_id)
        WHERE {f'u.acc="{accs}"'
            if isinstance(accs, str) else
            f"u.acc IN {accs}"}
        AND a.above_{threshold}=(1)
        GROUP BY inferred_compound_id, inferred_uniprot_acc
    '''

    # records = mysql_query(inferred_query, existing_conn=existing_conn)

    create_combined_sql = f'''
        SELECT 
            CASE
                WHEN predictions.predicted_compound_id IS NOT NULL THEN  predictions.predicted_compound_id
                ELSE predictions.inferred_compound_id
            END AS `compound_id`,
            CASE 
                WHEN predictions.predicted_uniprot_acc IS NOT NULL THEN  predictions.predicted_uniprot_acc
                ELSE predictions.inferred_uniprot_acc
            END AS `acc`,
            CASE 
                WHEN predictions.predicted_confidence_score IS NOT NULL AND 
                    predictions.inferred_confidence_score IS NOT NULL THEN 1000
                WHEN predictions.predicted_confidence_score IS NOT NULL
                    THEN predictions.predicted_confidence_score
                ELSE predictions.inferred_confidence_score
            END AS `confidence_score`,
            CASE 
                WHEN predictions.predicted_confidence_score IS NOT NULL AND 
                    predictions.inferred_confidence_score IS NOT NULL THEN "CONFIRMED"
                WHEN predictions.predicted_confidence_score IS NOT NULL
                    THEN "PREDICTED"
                ELSE "INFERRED"
            END AS `confidence_type`,
            predictions.all_targets
        FROM (
            WITH predicted AS ({predicted_query}), inferred AS ({inferred_query})
            SELECT *
            FROM predicted
            LEFT JOIN 
            inferred
                ON (predicted.predicted_uniprot_acc=inferred.inferred_uniprot_acc
                    AND predicted.predicted_compound_id=inferred.inferred_compound_id)
            UNION ALL
            SELECT * 
            FROM predicted
            RIGHT JOIN 
            inferred
                ON (predicted.predicted_uniprot_acc=inferred.inferred_uniprot_acc
                    AND predicted.predicted_compound_id=inferred.inferred_compound_id)
            WHERE predicted.predicted_uniprot_acc IS NULL
        ) AS predictions
      ORDER BY confidence_score DESC, confidence_type ASC
    '''

    # combined_table_name = "temp_compound"
    # create_temp_combined_table_sql = f'''
    # CREATE TEMPORARY TABLE {combined_table_name} (
    #     compound_id MEDIUMINT, 
    #     acc VARCHAR(10),
    #     confidence_score SMALLINT,
    #     confidence_type VARCHAR(10),
    #     all_targets VARCHAR(1000),
    #     PRIMARY KEY(compound_id)
    # )
    # {create_combined_sql}
    # '''
    # mysql_create_table(create_temp_combined_table_sql, existing_conn=existing_conn)

    query = f'''
        SELECT c.coconut_id AS `id`, 
        c.image AS `image`, 
        c.name AS `name`, 
        c.formula AS `formula`, 
        c.smiles AS `smiles`,
        predictions.acc,
        predictions.confidence_score,
        predictions.confidence_type,
        predictions.all_targets
        FROM compounds AS c
        INNER JOIN ({create_combined_sql}) AS `predictions`
            ON (c.compound_id=predictions.compound_id)
        ORDER BY confidence_score DESC, confidence_type ASC  
    '''

    records, cols = mysql_query(query, return_cols=True, existing_conn=existing_conn)
    existing_conn.close()

    if as_dict:
        records = [
            {c:r for c, r in zip(cols, record)}
            for record in records
        ]
    return records, cols

def get_combined_uniprot_confidences_for_compounds(
    compound_ids, 
    threshold=750, 
    as_dict=True,
    existing_conn=None):

    if isinstance(compound_ids, list) or isinstance(compound_ids, set):
        compound_ids = tuple(compound_ids)
    else:
        assert isinstance(compound_ids, str)

    if existing_conn is None:
        existing_conn = connect_to_mysqldb()

    predicted_query = f'''
        SELECT c.coconut_id AS `predicted_coconut_id`, 
            cu.uniprot_id AS `predicted_uniprot_id`, 
            cu.confidence_score AS `predicted_confidence_score`
        FROM compounds AS c
        INNER JOIN compound_to_uniprot AS cu 
            ON (c.compound_id=cu.compound_id)
        WHERE {f'c.coconut_id="{compound_ids}"'
            if isinstance(compound_ids, str) else
            f"c.coconut_id IN {compound_ids}"}
        AND cu.above_{threshold}=(1)
    '''

    inferred_query = f'''
        SELECT c.coconut_id AS `inferred_coconut_id`, 
            tu.uniprot_id AS `inferred_uniprot_id`, 
            MAX(a.confidence_score) AS `inferred_confidence_score`,
            JSON_ARRAYAGG(
                JSON_OBJECT(
                    'target_name', t.target_name,
                    'confidence_score', a.confidence_score
                )
            ) AS all_targets
        FROM compounds AS c
        INNER JOIN activities AS a
            ON (c.compound_id=a.compound_id)
        INNER JOIN targets AS t 
            ON (a.target_id=t.target_id)
        INNER JOIN targets_to_uniprot as tu 
            on (t.target_id=tu.target_id) 
        WHERE {f'c.coconut_id="{compound_ids}"'
            if isinstance(compound_ids, str) else
            f"c.coconut_id IN {compound_ids}"}
        AND a.above_{threshold}=(1)
        GROUP BY inferred_coconut_id, inferred_uniprot_id
    '''

    create_combined_sql = f'''
        SELECT 
            CASE
                WHEN predictions.predicted_coconut_id IS NOT NULL THEN  predictions.predicted_coconut_id
                ELSE predictions.inferred_coconut_id
            END AS `compound_id`,
            CASE 
                WHEN predictions.predicted_uniprot_id IS NOT NULL THEN  predictions.predicted_uniprot_id
                ELSE predictions.inferred_uniprot_id
            END AS `uniprot_id`,
            CASE 
                WHEN predictions.predicted_confidence_score IS NOT NULL AND 
                    predictions.inferred_confidence_score IS NOT NULL THEN 1000
                WHEN predictions.predicted_confidence_score IS NOT NULL
                    THEN predictions.predicted_confidence_score
                ELSE predictions.inferred_confidence_score
            END AS `confidence_score`,
            CASE 
                WHEN predictions.predicted_confidence_score IS NOT NULL AND 
                    predictions.inferred_confidence_score IS NOT NULL THEN "CONFIRMED"
                WHEN predictions.predicted_confidence_score IS NOT NULL
                    THEN "PREDICTED"
                ELSE "INFERRED"
            END AS `confidence_type`,
            predictions.all_targets AS `all_targets`
        FROM (
            WITH predicted AS ({predicted_query}), inferred AS ({inferred_query})
            SELECT *
            FROM predicted
            LEFT JOIN 
            inferred
                ON (predicted.predicted_uniprot_id=inferred.inferred_uniprot_id
                    AND predicted.predicted_coconut_id=inferred.inferred_coconut_id)
            UNION ALL
            SELECT * 
            FROM predicted
            RIGHT JOIN 
            inferred
                ON (predicted.predicted_uniprot_id=inferred.inferred_uniprot_id
                    AND predicted.predicted_coconut_id=inferred.inferred_coconut_id)
            WHERE predicted.predicted_uniprot_id IS NULL
        ) AS predictions
    '''
   
    # combined_table_name = "temp_uniprot"
    # create_temp_combined_table_sql = f'''
    # CREATE TEMPORARY TABLE {combined_table_name} (
    #     compound_id VARCHAR(10), 
    #     uniprot_id MEDIUMINT,
    #     confidence_score SMALLINT,
    #     confidence_type VARCHAR(10),
    #     all_targets VARCHAR(1000),
    #     PRIMARY KEY(uniprot_id)
    # )
    # {create_combined_sql}
    # '''
    # mysql_create_table(create_temp_combined_table_sql, existing_conn=existing_conn)

    query = f'''
        select predictions.compound_id,
        u.acc,
        u.protein,
        u.gene,
        u.organism_common,
        u.organism_scientific,
        u.organism_synonym,
        predictions.confidence_score,
        predictions.confidence_type,
        predictions.all_targets
        FROM ({create_combined_sql}) AS predictions
        INNER JOIN uniprot AS u 
            ON (predictions.uniprot_id=u.uniprot_id)
        ORDER BY confidence_score DESC, confidence_type ASC
    '''

    records, cols = mysql_query(query, return_cols=True, existing_conn=existing_conn)
    existing_conn.close()

    if as_dict:
        records = [
            {c:r for c, r in zip(cols, record)}
            for record in records
        ]
    return records, cols

if __name__ == "__main__":

    from natural_products.views import VALID_ORGANISMS

    # records, cols = get_info_for_multiple_compounds(columns=("compound_id", "kingdom_name", "species_name"), kingdom_name="Marine")
    # records = get_all_species(species_group="Aconitum")
    # records = get_all_reactions(organisms=VALID_ORGANISMS)

    # for record in records:
    #     if record[1] != "Homo sapiens":
    #         print (record)

    # print (len(records))

    # acc = "O42713"
    # records, cols = get_drugs_for_uniprots(acc)
    # records, cols = get_diseases_for_uniprots(acc)


    # # records, cols = get_diseases_for_drugs("N4-(3-chlorophenyl)quinazoline-4,6-diamine")

    # for record in records[:5]:
    #     print (record)

    # print (len(records))
    # # pd.DataFrame(records, columns=cols).to_csv("drugs.csv")

    pathways = "Signaling by EGFR"
    # reaction = "ERBB4 forms heterodimers with EGFR"
    organism = "Homo sapiens"

    reaction = [
     "VEGFA dimer:p-6Y-VEGFR2 dimer:PI3K phosphorylates PIP2 to PIP3",                      
    "p-6Y-VEGFR2 binds PI3K",                                                              
    "ERBB2 forms heterodimers with ligand-activated ERBB receptors: EGFR, ERBB3 and ERBB4",
    "ERBB4 forms heterodimers with EGFR",                                                   
    "NRP-1 forms a ternary complex with VEGF165 and VEGFR1",
    ]

    from timeit import default_timer

    start_time = default_timer()

    # records, cols = get_pathway_hits(pathways, organism, threshold=750, 
    #     min_target_coverage=0.1, min_pathway_coverage=0)
    records, cols = get_reaction_hits(reaction, organism, threshold=950, min_target_coverage=.1, min_reaction_coverage=.2)

    # print (records[0])
    print (len(records))

    print (default_timer() - start_time)

    # coconut_id = "CNP0000002"
    # threshold = 750

    # # records, cols = get_combined_uniprot_confidences_for_compounds(coconut_id, threshold=threshold, as_dict=True)
    # records, cols = get_combined_compound_confidences_for_uniprots("P08183", threshold=threshold, as_dict=True)

    # for record in records[:5]:
    #     print (record)

    # print (len(records))

    # records = get_inferred_uniprots_for_compounds("CNP0000002")

    # accs = get_all_uniprot()

    # print(accs[:5])

    # acc = "P08183"
    # # records = get_combined_compound_confidences_for_uniprots(acc, threshold=950)
    # records = get_combined_uniprot_confidences_for_compounds("CNP0000002")
    # # records = get_inferred_compounds_for_uniprots(acc)
    # print (len(records))
    # records = get_combined_uniprot_confidences_for_compounds("CNP0000002", threshold=750)
    # combined_attempt_accs = {record[1] for record in records}
    # records = get_inferred_compounds_for_uniprots(acc)

    # with open(f"{acc}_hits.tsv", "w") as f:
    # for record in records:
        # print (record)
    #         f.write("\t".join(map(str, record)) + "\n")
    # num_combined = len(records)

    # records = get_predicted_uniprots_for_compounds("CNP0000002", threshold=750)
    # predicted_accs = {record[1] for record in records}

    # records = get_inferred_uniprots_for_compounds("CNP0000002", threshold=750)
    # inferred_accs = {record[1] for record in records}

    # combined_accs = predicted_accs.union(inferred_accs)
    # intersect_accs = predicted_accs.intersection(inferred_accs)
  
    # print (sorted(predicted_accs)[:5])
    # print (sorted(inferred_accs)[:5])
    # print (sorted(combined_accs)[:5])

    # print (num_combined)

    # print (len(predicted_accs))
    # print (len(inferred_accs))
    # print (len(combined_accs))
    # print (len(intersect_accs))
    # print (intersect_accs)

    # print (len(predicted_accs.intersection(combined_attempt_accs)))
    # print (len(inferred_accs.intersection(combined_attempt_accs)))

    # hits = get_all_activities_for_compound("CNP0000002", threshold=650)

    # print (hits["All_targets"])

    # smiles = get_info_for_multiple_compounds(compounds=None,
    #     columns=("compound_id", "smiles"))

    # from utils.io import write_smiles
    # write_smiles(smiles, "coconut_smiles.smi")

    # import pandas as pd
    # targets = pd.read_csv("models/target_ids.txt", 
    #     index_col=0, header=None, names=["uniprot"])
    # targets = targets["uniprot"]

    # query = f'''
    #     SELECT acc, uniprot_id
    #     FROM uniprot
    #     WHERE acc IN {tuple(targets)} 
    # '''

    # records = mysql_query(query)

    # acc_to_db_id = {acc: _id for acc, _id in records}

    # id_to_db_id = {_id: acc_to_db_id[acc]
    #     for _id, acc in targets.items()}

    # import json 
    # with open("id_to_db_id.json", "w") as f:
    #     json.dump(id_to_db_id, f, sort_keys=True, indent=4)
   
    # compound_id=["CNP0000002", "CNP0000005"]

    # records = get_uniprots_for_compound(compound_id, threshold=700)
    # for record in records:
    #     print (record)


    # name_like = "gink"

    # records = get_info_for_multiple_compounds(name_like=None, compound_ids=None, formula_like="H%2")

    # for record in records:
    #     print (record)