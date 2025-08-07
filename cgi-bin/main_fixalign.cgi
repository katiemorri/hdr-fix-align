#!/usr/bin/env python3

#FixAlign - Analyze sequence alignments, manually approve mismatches

import cgi
import jinja2
from Bio import SeqIO
from Bio.Align import PairwiseAligner
import cgitb
cgitb.enable()

#print header, generate html in browser
print("Content-Type: text/html\n\n")

#show jinja2 template to the loader (where template file is located);create env and load template
template_loader = jinja2.FileSystemLoader("/var/www/html/kmorri74/fixalign_proj/templates")
env = jinja2.Environment(loader=template_loader)
template = env.get_template("results_display.html")

#access from data (use placeholders)
form = cgi.FieldStorage()

#parse input html form
#make sure no None return for string, make sequences homogenous (no white space inside or outside string)
subject_seq = form.getvalue("subject_seq", "").replace(" ","").replace("\n", "").strip().upper()
query_seq = form.getvalue("query_seq", "").replace(" ", "").replace("\n", "").strip().upper()
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

#perform alignment via pairwisealigner(new-not pairwise2), biopython
#gap penalties are preset in this script to -4 and -1 (open, and then penalty)
aligner = PairwiseAligner()
aligner.mode ='global'
aligner.match_score = match_score
aligner.mismatch_score = mismatch_score
aligner.open_gap_score = -4
aligner.extend_gap_score = -1
alignments = aligner.align(subject_seq, query_seq) #perform alignment

#assign alignment variables
initial_alignment = alignments[0]
alignment_string = initial_alignment.format().splitlines()
aligned_subject = alignment_string[0]
final_matched_lines = alignment_string[1]
aligned_query = alignment_string[2]

#aligned_subject = initial_alignment.target
#aligned_query = initial_alignment.query
score = initial_alignment.score

coordinates = initial_alignment.coordinates
subject_start = coordinates[0,0]
subject_end = coordinates[0,-1]
query_start = coordinates[1,0]
query_end = coordinates[1,-1]

#identify mismatch & gap positions (misalignments)
valid_count_chars = ['|', '-', '.']
start_index = 0
for i, char in enumerate(final_matched_lines):
    if char in valid_count_chars:
        start_index = i
        break

total_bases = len(final_matched_lines) - start_index
base_matches = 0
misaligned_positions = []

#wrap mismatch positions for highlighting (CSS) and checkbox in JS
highlighted_subject = []
highlighted_query = []

for i in range(len(final_matched_lines)):
    subject_base = aligned_subject[i]
    query_base = aligned_query[i]
    match_char = final_matched_lines[i]

    #check if not aligned, then save position
    if match_char == '|':               #match
        base_matches += 1               #save for calc. percentage
        #wappend matches without wrap
        highlighted_subject.append(subject_base)
        highlighted_query.append(query_base)
    else:
        misaligned_positions.append(i)  #mis-match
        #append mismatch to have the wrap
        highlighted_subject.append(f"<span class='mismatch'>{subject_base}</span>")
        highlighted_query.append(f"<span class='mismatch'>{query_base}</span>")

#calc percentage of alignment between subject and query (including gaps)
if total_bases > 0:
    align_percent = round((base_matches / total_bases) * 100, 2)
else:
    align_percent = 'N/A'

#total # of misalignments
total_misalignments = total_bases - base_matches

#create function to block aligned highlighted sub, match lines, and highlighted, aligned seq
def blocked_lines(seq, block=60):
    return [''.join(seq[i:i+block]) for i in range(0, len(seq), block)]
blocked_matched_lines = blocked_lines(final_matched_lines)
highlighted_subject_lines = blocked_lines(highlighted_subject)
highlighted_query_lines = blocked_lines(highlighted_query)

#determine pass/fail based on % match thresholds
pass_status = "Fail"
if clone == "Forward":
    if arm_type == "5Arm":
        if align_percent >= 95:
            pass_status = "Pass"
    elif arm_type == "3Arm":
        if align_percent >= 98:
            pass_status = "Pass"
elif clone == "Reverse":
    if arm_type == "3Arm":
        if align_percent >= 95:
            pass_status = "Pass"
    elif arm_type == "5Arm":
        if align_percent >= 98:
            pass_status = "Pass"
else:
    pass_status = "Fail"

print(template.render(
    subject_seq=subject_seq, #from html input form
    query_seq=query_seq,
    scoring=scoring,
    arm_type=arm_type,
    clone=clone,
    aligned_subject=aligned_subject, #from main cgi
    aligned_query=aligned_query,
    score=score,
    subject_start=subject_start,
    subject_end=subject_end,
    query_start=query_start,
    query_end=query_end,
    align_percent=align_percent,
    misaligned_positions=misaligned_positions,
    total_misalignments=total_misalignments,
    final_matched_lines=final_matched_lines,
    pass_status=pass_status,
    highlighted_subject_lines=highlighted_subject_lines,
    blocked_matched_lines=blocked_matched_lines,
    highlighted_query_lines=highlighted_query_lines))
