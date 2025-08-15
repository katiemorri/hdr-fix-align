#!/usr/bin/env python3
import cgi
import mysql.connector
import jinja2
import cgitb
cgitb.enable()

#print header, generate html in browser
print("Content-Type: text/html\n\n")

#show jinja2 template to the loader (where template file is located);create env and load template
template_loader = jinja2.FileSystemLoader("/var/www/html/kmorri74/fixalign_proj/templates")
env = jinja2.Environment(loader=template_loader)
template = env.get_template("database_confirm.html")


# Get form data
form = cgi.FieldStorage()
subject_seq = form.getvalue('subject_seq')
query_seq = form.getvalue('query_seq')
pairwise_score = form.getvalue('pairwise_score')
scoring = form.getvalue('scoring')
alignment_score = form.getvalue('alignment_score')
pass_status = form.getvalue('pass_status')

# Remove % to store
if alignment_score and alignment_score.endswith('%'):
    alignment_score = alignment_score.replace('%','').strip()

# Connect to MySQL
conn = mysql.connector.connect(
    host="localhost",
    user="kmorri74",
    password="MercuryMob$001",
    database="kmorri74_chado"
)
cursor = conn.cursor()

# Insert into the database with cursor execute
sql = '''
INSERT INTO Alignment_Results (
    subject_seq,
    query_seq,
    pairwise_score,
    scoring,
    alignment_score,
    pass_status)
    VALUES (%s, %s, %s, %s, %s, %s)
'''
cursor.execute(sql, (subject_seq, query_seq, pairwise_score, scoring, alignment_score, pass_status))
conn.commit()
cursor.close()
conn.close()

print(template.render(
    subject_seq=subject_seq,
    query_seq=query_seq,
    pairwise_score=pairwise_score,
    scoring=scoring,
    alignment_score=alignment_score,
    pass_status=pass_status,
    confirmation="Results were successfully saved to database."
))
