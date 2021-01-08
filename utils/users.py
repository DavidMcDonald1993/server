import sys
import os.path
sys.path.insert(1, 
    os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))

import secrets

from utils.mysql_utils import mysql_create_table, mysql_query, mysql_insert_many
from utils.email_utils import send_mail

def get_token(token_length=128): # token is double length
    return secrets.token_hex(token_length)

def create_files_table():

    create_files_table_sql = '''
        CREATE TABLE user_files 
        (
            id int NOT NULL UNIQUE AUTO_INCREMENT,
            token VARCHAR(300) UNIQUE NOT NULL, 
            user_id INT NOT NULL, 
            path VARCHAR(255) NOT NULL,
            PRIMARY KEY (id),
            FOREIGN KEY(user_id) REFERENCES auth_user(id)
        )
        
    '''

    mysql_create_table(create_files_table_sql)

    return 0

def add_file_to_database(user_id, path, existing_conn=None):
    token = get_token()

    print (len(token))

    insert_sql = '''
        INSERT INTO user_files(token, user_id, path)
        VALUES (%s, %s, %s)
    '''
    
    rows = [(token, user_id, path)]

    mysql_insert_many(insert_sql, rows, existing_conn=existing_conn)

    return token

def get_file_from_token(token, user_id, existing_conn=None):

    query = f'''
        SELECT path
        FROM user_files
        WHERE token="{token}"
        AND user_id={user_id}
        LIMIT 1
    '''
    
    records = mysql_query(query, existing_conn=existing_conn)
    if len(records) == 1:
        return records[0][0]
    else:
        return None

def send_file_to_user(user, filename, max_MB=5):
    if os.path.getsize(filename) / (1024*1024) < max_MB: # file smaller than max MB
        # send as attachment
        send_mail(user.username,
            user.email,
            attach_file_name=filename)
    else:
        
        # file is too large: save on server for download later
        token = add_file_to_database(user.id, path=filename)
        print ("added file", filename, "to database for user", user.username)
        print ("generated token", token)
        send_mail(user.username,
            user.email,
            token=token)

    # TODO delele file?
    # os.remove(filename)

    return 0

if __name__ == "__main__":

    # token = add_file_to_database(1, path="test.sdf")

    token = "f5cfd4845379b498755a7fa4f683776be8c1aa8949aa219809f78bec7c0ee59532561ab033421417fea30d81437538ec5fe47cfd4881e4dd16009e0676812a954855210daeb13f0f0f9f6136cbaa82814dd63b2cb8998e41b8fcc434eeb64446d208e0146fa94393fd7cbdf8c406a04a6fca069ae0fa668caaeade3b35197311"

    print (token)

    path = get_file_from_token(token, user_id=1)

    print (path)