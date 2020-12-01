import os

import smtplib
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from email.mime.base import MIMEBase
from email import encoders

import json 

import argparse

with open("email_credentials.json", "r") as f:
    credentials = json.load(f)

def send_mail(
    receiver_name,
    receiver_address,
    attach_file_name,
    subject="NPAIEngine results",
    mail_content="Dear {},\n\nPlease find attached the results from NPAIEngine.\n\nSincerely,\nShan He"
    ):

    print ("sending email to user", receiver_name, "at address", receiver_address)
    print ("sending file", attach_file_name)

    #Setup the MIME
    message = MIMEMultipart()
    message['From'] = credentials["sender_address"]
    message['To'] = receiver_address
    message['Subject'] = subject

    #The body and the attachments for the mail
    mail_content = mail_content.format(receiver_name)
    message.attach(MIMEText(mail_content, 'plain'))
    attach_file = open(attach_file_name, 'rb') # Open the file as binary mode
    payload = MIMEBase('application', 'octate-stream')
    payload.set_payload((attach_file).read())
    encoders.encode_base64(payload) #encode the attachment
 
    #add payload header with filename
    payload.add_header("Content-disposition", 'attachment', filename=os.path.basename(attach_file_name))
    message.attach(payload)
    
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
    "/home/david/Desktop/test-PASS-out.sdf")