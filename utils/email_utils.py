import os

import sys
import os.path
sys.path.insert(1, 
    os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))

import smtplib
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from email.mime.base import MIMEBase
from email import encoders

import json 

import argparse

from utils.io import load_json

def send_mail(
    receiver_name,
    receiver_address,
    attach_file_name=None,
    token=None,
    subject="NPAIEngine results",
    # mail_content="Dear {},\n\nPlease find attached the results from NPAIEngine.\n\nSincerely,\nShan He"
    ):

    host_information = load_json("host_credentials.json")
    credentials = load_json("email_credentials.json")

    print ("sending email to user", receiver_name, "at address", receiver_address, "with subject", subject)
    print ("sending file", attach_file_name)

    #Setup the MIME
    message = MIMEMultipart()
    message['From'] = credentials["sender_address"]
    message['To'] = receiver_address
    message['Subject'] = subject

    #The body and the attachments for the mail
    if attach_file_name is not None:
        mail_content=f'''
        Dear {receiver_name.capitalize()},
        
        Please find attached the results from NPAIEngine.
        
        Sincerely,
        Shan He
        '''

        message.attach(MIMEText(mail_content, 'plain'))
        attach_file = open(attach_file_name, 'rb') # Open the file as binary mode
        payload = MIMEBase('application', 'octate-stream')
        payload.set_payload((attach_file).read())
        encoders.encode_base64(payload) #encode the attachment
    
        #add payload header with filename
        payload.add_header("Content-disposition", 'attachment', filename=os.path.basename(attach_file_name))
        message.attach(payload)
    else:
        mail_content=f'''
        Dear {receiver_name},
        
        Please find your files at http://{host_information["host"]}:{host_information["port"]}/download/{token}
        
        Sincerely,
        Shan He
        '''

        message.attach(MIMEText(mail_content, 'plain'))

    
    #Create SMTP session for sending the mail
    session = smtplib.SMTP('smtp.gmail.com', 587) #use gmail with port
    session.starttls() #enable security
    session.login(credentials["sender_address"], credentials["sender_pass"]) #login with mail_id and password
    text = message.as_string()
    session.sendmail(credentials["sender_address"], receiver_address, text)
    session.quit()
    print('Mail sent')

if __name__ == "__main__":
    send_mail("David", "davemcdonald93@gmail.com", 
    # attach_file_name="/home/david/Desktop/test-PASS-out.sdf"
    token="f5cfd4845379b498755a7fa4f683776be8c1aa8949aa219809f78bec7c0ee59532561ab033421417fea30d81437538ec5fe47cfd4881e4dd16009e0676812a954855210daeb13f0f0f9f6136cbaa82814dd63b2cb8998e41b8fcc434eeb64446d208e0146fa94393fd7cbdf8c406a04a6fca069ae0fa668caaeade3b35197311"
    )