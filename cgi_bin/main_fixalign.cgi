#!/user/bin/env python

# FixAlign - Analyze sequence aligments, and easily approve mismatches

import cgi
import jinja2
from Bio import SeqIO
from Bio.pairwise2 import format_alignment

#print header, generate html in browser
print("Content-Type: text/html\n\n")

#show jinja2 template to the loader (where template file is located);create env and load template
template_loader = jinja2.FileSystemLoader(/'')
env = jinja2.Environment(loader=template_loader)
template = env.get_template("results_display.html")

#access from data (use placeholders)
form = cgi.FieldStorage()

#parse input html form
#make sure no None return for string, make sequences homogenous (no white space inside or outside string)
subject_seq = form.getvalue("subject_seq", "").replace(" ","").replace("\n", "").strip().upper
query_seq = form.getvalue("query_seq", ""),replace(" ", "").replace("\n", "").strip().upper
scoring = form.getvalue("scoring", "1,0").strip()
arm_type = form.getvalue("arm", "").strip()
clone = form.getvalue("clone", "").strip()

    #parse scoring sring into match and mismatch score via coma; make intergers
    try:
        split_score = scoring.split(",")
        match_score = int(split_score[0].strip())
        mismatch_score = int(split_score[1].strip())
    #error of most sort, will be default of (1,0)
    except:
        match_score = 1
        mismatch_score = 0

#perform alignment via pairwise2, biopython
#gap penalties are preset in this script to -4 and -1 (open, and then penalty)
alignments = pairwise2.align.globalms(subject_seq, query_seq, match_score, mismatch_score, -4, -1)
#use first alignment made
initial_aligment = alignments[0]

#visual alignment

#initial (unedited) alignment score

#frameshift detection (Y/N)

#total # of  mismatches

#5 Arm vs 3 Arm region analysis (if more misalignments occur in the 3 Arm vs the 5 Arm region

#Determine pass/fail based on % match thresholds


#print(template.render(
