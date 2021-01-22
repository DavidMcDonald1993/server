import sys
import os.path
sys.path.insert(1, 
    os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))

from collections import defaultdict

from utils.mysql_utils import mysql_query, sanitise_names

import urllib.parse as urlparse

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
        if not isinstance(targets, tuple):
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
    SELECT c.coconut_id AS `ID`, cu.confidence_score AS `Confidence Score`, 
        u.acc AS `Uniprot ACC`
    FROM compounds AS c
    INNER JOIN compound_to_uniprot AS cu 
        ON (c.compound_id=cu.compound_id)
    INNER JOIN uniprot AS u 
        ON (cu.uniprot_id=u.uniprot_id)
    WHERE cu.above_{threshold}=(1)
    AND {f"c.coconut_id='{coconut_id}'" if isinstance(coconut_id, str)
        else f"c.coconut_id IN {coconut_id}"}
    '''

    # get_uniprots_sql = f'''
    #     SELECT c.coconut_id, t.target_name, u.acc
    #     FROM compounds AS c
    #     INNER JOIN activities AS a ON (c.compound_id=a.compound_id)
    #     INNER JOIN targets AS t ON (a.target_id=t.target_id)
    #     INNER JOIN targets_to_uniprot AS tu ON (a.target_id=tu.target_id)
    #     INNER JOIN uniprot AS u ON (tu.uniprot_id=u.uniprot_id)
    #     WHERE {f"c.coconut_id='{coconut_id}'" if isinstance(coconut_id, str)
    #         else f"c.coconut_id IN {coconut_id}"}
    #     {f"AND a.Pa>{threshold}" if threshold>0 else ""}
    #     {f"AND a.Pa>a.Pi" if filter_pa_pi else ""}
    # '''

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
    SELECT GROUP_CONCAT(DISTINCT(u.acc) SEPARATOR '-') AS `Uniprot ACCS`, 
        COUNT(DISTINCT(u.acc)) AS `Number Uniprot ACCs`,
        p.num_uniprot AS `Total Uniprot ACCs`,
        COUNT(DISTINCT(u.acc)) / p.num_uniprot AS `Uniprot Coverage`,      
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
    records = mysql_query(query, existing_conn=existing_conn)

    records = [
        (uniprots, uniprot_count, total_uniprots,
            coverage, pathway_name, urlparse.quote(pathway_name),
            organism, urlparse.quote(organism), url)
        for uniprots, uniprot_count, total_uniprots,
            coverage, pathway_name,
            organism, url in records
    ]

    return records

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
    SELECT GROUP_CONCAT(DISTINCT(u.acc) SEPARATOR '-') AS `Uniprot ACCS`, 
        COUNT(DISTINCT(u.acc)) AS `Number Uniprot ACCs`,
        r.num_uniprot AS `Total Uniprot ACCs`,
        COUNT(DISTINCT(u.acc)) / r.num_uniprot AS `Uniprot Coverage`,   
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
    records = mysql_query(query, existing_conn=existing_conn)

    records = [
        (uniprots, uniprot_count, total_uniprots,
            coverage, reaction_name, urlparse.quote(reaction_name,),
            organism, urlparse.quote(organism), url)
        for uniprots, uniprot_count, total_uniprots,
            coverage, reaction_name,
            organism, url in records
    ]

    return records


def get_target_hits(
    targets, 
    thresholds,
    filter_pa_pi=True,
    limit=None
    ):
    assert filter_pa_pi
    assert isinstance(targets, list)
    if not isinstance(thresholds, list):
        assert isinstance(thresholds, int), thresholds
        # assert thresholds in (900, )
        thresholds = [thresholds]
    num_targets = len(targets)
    if len(thresholds) < num_targets:
        thresholds = thresholds * num_targets

    print ("querying PASS database for compounds",
        "that hit targets", targets,
         "with Pa greater than or equal to",
        "thresholds", thresholds)

    target_names = sanitise_names(targets)
    columns = ", ".join((
        f'''
            `{target}_activity`.Pa AS `{target}-Pa`, `{target}_activity`.Pi AS `{target}-Pi`, 
            `{target}_activity`.confidence_score AS `{target}-Confidence Score`
        '''
        for target in target_names
    ))
    tables = "\n".join((
            f'''
            INNER JOIN activities AS `{target}_activity` 
                ON (c.compound_id=`{target}_activity`.compound_id)
            INNER JOIN targets AS `{target}_target` 
                ON (`{target}_activity`.target_id=`{target}_target`.target_id)
            '''
        for target in target_names
    ))
    conditions = "\n".join((
        f'''
        AND `{target_name}_activity`.above_{threshold}=(1) 
        AND `{target_name}_target`.target_name="{target}"
        '''
        for target_name, target, threshold in zip(target_names[1:], targets[1:], thresholds[1:])      
    ))
    query = f'''
        SELECT c.compound_id, c.coconut_id AS `ID`, c.image_path AS `Image`, c.name AS `Molecule Name`, 
            c.formula AS `Molecular Formula`, c.clean_smiles AS `SMILES`,
        {columns}
        FROM compounds AS c
        {tables}
        WHERE `{target_names[0]}_activity`.above_{thresholds[0]}=(1)
        AND `{target_names[0]}_target`.target_name="{targets[0]}"
        {conditions}
        ORDER BY `{target_names[0]}-Confidence Score` DESC 
        {f"LIMIT {limit}" if limit is not None else ""}
    '''

    hits, cols = mysql_query(query, return_cols=True)

    records = [
        record[1:] for record in hits
    ]

    return records, cols[1:]

def get_pathway_hits(
    pathways,
    organisms, 
    threshold=0, 
    filter_pa_pi=True, 
    limit=None,
    existing_conn=None):
    assert filter_pa_pi

    if not isinstance(pathways, list):
        assert isinstance(pathways, str)
        pathways = [pathways]
    if not isinstance(organisms, list):
        assert isinstance(organisms, str)
        organisms = [organisms]

    pathway_names = sanitise_names(
        [   f"{pathway}_{organism}"
            for pathway, organism in zip(pathways, organisms)]
    )

    columns = ", ".join((
        f'''
        GROUP_CONCAT(DISTINCT(`{pathway}_target`.target_name) SEPARATOR '-') AS `{pathway} Target Names`, 
        COUNT(DISTINCT(`{pathway}_target`.target_name)) AS `{pathway} Number Targets`, 
        GROUP_CONCAT(DISTINCT(`{pathway}_uniprot`.acc) SEPARATOR '-') AS `{pathway} Uniprot ACCs`, 
        COUNT(DISTINCT(`{pathway}_uniprot`.acc)) AS `{pathway} Number Uniprot ACCs`, 
        `{pathway}_pathway`.num_uniprot AS `{pathway} Total Uniprot ACCs`,
        CAST(COUNT(DISTINCT(`{pathway}_uniprot`.acc)) / `{pathway}_pathway`.num_uniprot AS CHAR)
            AS `{pathway} Uniprot Coverage`,   
        `{pathway}_pathway`.pathway_name AS `{pathway} Name`,
        `{pathway}_pathway`.organism AS `{pathway} Organism`,
        `{pathway}_pathway`.pathway_url AS `{pathway} URL`
        '''
        for pathway in pathway_names
    ))

    tables = "\n".join((
        f'''
        INNER JOIN activities AS `{pathway}_activity` 
            ON (c.compound_id=`{pathway}_activity`.compound_id)
        INNER JOIN targets AS `{pathway}_target` 
            ON (`{pathway}_activity`.target_id=`{pathway}_target`.target_id)
        INNER JOIN targets_to_uniprot AS `{pathway}_targets_to_uniprot` 
            ON (`{pathway}_activity`.target_id=`{pathway}_targets_to_uniprot`.target_id)
        INNER JOIN uniprot AS `{pathway}_uniprot` 
            ON (`{pathway}_targets_to_uniprot`.uniprot_id=`{pathway}_uniprot`.uniprot_id)
        INNER JOIN uniprot_to_pathway AS `{pathway}_uniprot_to_pathway` 
            ON (`{pathway}_targets_to_uniprot`.uniprot_id=`{pathway}_uniprot_to_pathway`.uniprot_id)
        INNER JOIN pathway AS `{pathway}_pathway`
            ON (`{pathway}_uniprot_to_pathway`.pathway_id={pathway}_pathway.pathway_id)
        '''
        for pathway in pathway_names
    ))

    conditions = "\n".join((
        f'''
        AND `{pathway_name}_activity`.above_{threshold}=(1) 
        AND `{pathway_name}_pathway`.pathway_name="{pathway}"
        AND `{pathway_name}_pathway`.organism="{organism}"
        '''
        for pathway_name, pathway, organism in zip(pathway_names[1:], pathways[1:], organisms[1:])
    ))

    group_by = ",".join((
        f'''
            `{pathway_name}_pathway`.pathway_name,
            `{pathway_name}_pathway`.organism,
            `{pathway_name}_pathway`.pathway_url
       '''
        for pathway_name in pathway_names
    ))

    query = f'''
        SELECT c.compound_id, c.coconut_id AS `ID`, c.image_path AS `Image`, c.name AS `Molecule Name`, 
            c.formula AS `Molecular Formula`, c.clean_smiles AS `SMILES`,
        {columns}
        FROM compounds AS c
        {tables}
        WHERE `{pathway_names[0]}_activity`.above_{threshold}=(1)
        AND `{pathway_names[0]}_pathway`.pathway_name="{pathways[0]}"
        AND `{pathway_names[0]}_pathway`.organism="{organisms[0]}"
        {conditions}
        GROUP BY c.compound_id, c.coconut_id, c.name, c.formula, c.clean_smiles, {group_by}
        ORDER BY `{pathway_names[0]} Uniprot Coverage` DESC
        {f"LIMIT {limit}" if limit is not None else ""}
    '''

    records, cols = mysql_query(query, return_cols=True,)

    records = [
        record[1:] for record in records
    ]
    return records, cols[1:]

def get_reaction_hits(
    reactions,
    organisms,
    threshold=0,
    filter_pa_pi=True, 
    limit=None,
    existing_conn=None):
    assert filter_pa_pi

    if not isinstance(reactions, list):
        assert isinstance(reactions, str)
        reactions = [reactions]
    if not isinstance(organisms, list):
        assert isinstance(organisms, str)
        organisms = [organisms]

    reaction_names = sanitise_names(
        [   f"{reaction}_{organism}"
            for reaction, organism in zip(reactions, organisms)]
    )
    columns = ", ".join((
        f'''
        GROUP_CONCAT(DISTINCT(`{reaction}_target`.target_name) SEPARATOR '-') AS `{reaction} Target Names`, 
        COUNT(DISTINCT(`{reaction}_target`.target_name)) AS `{reaction} Number Targets`, 
        GROUP_CONCAT(DISTINCT(`{reaction}_uniprot`.acc) SEPARATOR '-') AS `{reaction} Uniprot ACCs`, 
        COUNT(DISTINCT(`{reaction}_uniprot`.acc)) AS `{reaction} Number Uniprot ACCs`,
        `{reaction}_reaction`.num_uniprot AS `{reaction} Total Uniprot ACCs`,
        CAST(COUNT(DISTINCT(`{reaction}_uniprot`.acc)) / `{reaction}_reaction`.num_uniprot AS CHAR)
            AS `{reaction} Uniprot Coverage`,    
        `{reaction}_reaction`.reaction_name AS `{reaction} Name`,
        `{reaction}_reaction`.organism AS `{reaction} Organism`,
        `{reaction}_reaction`.reaction_url AS `{reaction} URL`
        '''
        for reaction in reaction_names
    ))

    tables = "\n".join((
        f'''
        INNER JOIN activities AS `{reaction}_activity` 
            ON (c.compound_id=`{reaction}_activity`.compound_id)
        INNER JOIN targets AS `{reaction}_target` 
            ON (`{reaction}_activity`.target_id=`{reaction}_target`.target_id)
        INNER JOIN targets_to_uniprot AS `{reaction}_targets_to_uniprot` 
            ON (`{reaction}_activity`.target_id=`{reaction}_targets_to_uniprot`.target_id)
        INNER JOIN uniprot AS `{reaction}_uniprot` 
            ON (`{reaction}_targets_to_uniprot`.uniprot_id=`{reaction}_uniprot`.uniprot_id)
        INNER JOIN uniprot_to_reaction AS `{reaction}_uniprot_to_reaction` 
            ON (`{reaction}_targets_to_uniprot`.uniprot_id=`{reaction}_uniprot_to_reaction`.uniprot_id)
        INNER JOIN reaction AS `{reaction}_reaction`
            ON (`{reaction}_uniprot_to_reaction`.reaction_id={reaction}_reaction.reaction_id)
        '''
        for reaction in reaction_names
    ))

    conditions = "\n".join((
        f'''
        AND `{reaction_name}_activity`.above_{threshold}=(1) 
        AND `{reaction_name}_reaction`.reaction_name="{reaction}"
        AND `{reaction_name}_reaction`.organism="{organism}"
        '''
        for reaction_name, reaction, organism in zip(reaction_names[1:], reactions[1:], organisms[1:])
    ))

    group_by = ",".join((
        f'''
        `{reaction_name}_reaction`.reaction_name,
        `{reaction_name}_reaction`.organism,
        `{reaction_name}_reaction`.reaction_url
       '''
        for reaction_name in reaction_names
    ))

    query = f'''
        SELECT c.compound_id, c.coconut_id AS `ID`, c.image_path AS `Image`, c.name AS `Molecule Name`, 
            c.formula AS `Molecular Formula`, c.clean_smiles AS `SMILES`,
        {columns}
        FROM compounds AS c
        {tables}
        WHERE `{reaction_names[0]}_activity`.above_{threshold}=(1)
        AND `{reaction_names[0]}_reaction`.reaction_name="{reactions[0]}"
        AND `{reaction_names[0]}_reaction`.organism="{organisms[0]}"
        {conditions}
        GROUP BY c.compound_id, c.coconut_id, c.name, c.formula, c.clean_smiles, {group_by}
        ORDER BY `{reaction_names[0]} Uniprot Coverage` DESC
        {f"LIMIT {limit}" if limit is not None else ""}
    '''

    records, cols = mysql_query(query, return_cols=True, )

    records = [
        record[1:] for record in records
    ]

    return records, cols[1:]

def get_info_for_multiple_compounds(
    compounds=None, 
    columns=("coconut_id", "name", "formula", "smiles"),
    limit=None):

    if compounds is not None:
        if not isinstance(compounds, str):
            if isinstance(compounds, list) or isinstance(compounds, set):
                compounds = tuple(compounds)
            assert isinstance(compounds, tuple)
            if len(compounds) == 1:
                compounds = compounds[0]

    compound_query = f'''
        SELECT {(", ".join(columns))}
        FROM compounds
        {f"WHERE coconut_id IN {compounds}" if isinstance(compounds, tuple)
        else f'WHERE coconut_id="{compounds}"' if isinstance(compounds, str) else ""}
        {f"LIMIT {limit}" if limit is not None else ""}
    '''

    return mysql_query(compound_query)

def get_categories():
    query = '''
    SELECT category_name
    FROM categories
    ORDER BY category_name ASC
    '''
    return [record[0] for record in mysql_query(query)]

def get_all_activities_for_compound(
    coconut_id, 
    threshold=0,
    filter_pa_pi=True):
    print ("getting all activities for compound", coconut_id)

    all_targets_query = f'''
    SELECT cat.category_name, t.target_name, a.Pa, a.Pi, a.confidence_score
    FROM categories AS cat
    INNER JOIN category_members AS cm ON (cat.category_id=cm.category_id)
    INNER JOIN targets AS t ON (cm.target_id=t.target_id) 
    INNER JOIN activities AS a ON (t.target_id=a.target_id)
    INNER JOIN compounds AS c ON (a.compound_id=c.compound_id)
    WHERE c.coconut_id='{coconut_id}'
    {f"AND a.above_{threshold}=(1)" if threshold > 0 and filter_pa_pi else ""}
    '''
    compound_hits = mysql_query(all_targets_query)

    compound_hits = [ # encode target name
        (category_name, target_name, urlparse.quote(target_name),
            pa, pi, confidence_score)
            for category_name, target_name, pa, pi, confidence_score
            in compound_hits
    ]
    
    category_activities = defaultdict(set)
    for category, *hit in compound_hits:
        hit = tuple(hit)
        category_activities["ALL_TARGETS"].add(hit)
        category_activities[category].add(hit)
   
    return dict(category_activities)

if __name__ == "__main__":

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
   
    compound_id=["CNP0000002", "CNP0000005"]

    records = get_uniprots_for_compound(compound_id, threshold=700)
    for record in records:
        print (record)
