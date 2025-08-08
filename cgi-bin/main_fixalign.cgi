#!/usr/bin/env python3

#FixAlign - Analyze sequence alignments, manually approve mismatches

import re
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

aligned_subject = ''
final_matched_lines = ''
aligned_query = ''

#put togther all lines in this block
triple_blocks = []
for i in range(0, len(alignment_string), 4):
    if i + 2 < len(alignment_string):
        subject_lines = alignment_string[i]
        matched_lines = alignment_string[i + 1]
        query_lines = alignment_string[i + 2]
        
#        triple_blocks.append(f"{subject_lines}\n{matched_lines}\n{query_lines}")
        
        highlighted_subject = []
        highlighted_match = []
        highlighted_query = []
        checkbox_lines = []

        for sub, mat, que in zip(subject_lines, matched_lines, query_lines):
            if mat == '|':
                highlighted_subject.append(f"<span>{sub}</span>")
                highlighted_match.append(f"<span>{mat}</span>")
                highlighted_query.append(f"<span>{que}</span>")
                checkbox_lines.append(' ')
            elif mat == '.' or mat == '-':
                highlighted_subject.append(f"<span class='mismatch'>{sub}</span>")
                highlighted_match.append(f"<span class='mismatch'>{mat}</span>")
                highlighted_query.append(f"<span class='mismatch'>{que}</span>")
                checkbox_lines.append("<input type='checkbox'>")
            else:
                highlighted_subject.append(sub)
                highlighted_match.append(mat)
                highlighted_query.append(que)
                checkbox_lines.append(' ')

        triple_blocks.append(
            ''.join(checkbox_lines) + "\n" +
            ''.join(highlighted_subject) + "\n" +
            ''.join(highlighted_match) + "\n" +
            ''.join(highlighted_query)
        )

        aligned_subject += subject_lines
        final_matched_lines += matched_lines
        aligned_query += query_lines

formatted_alignment = "\n\n".join(triple_blocks)

#aligned_subject = initial_alignment.target
#aligned_query = initial_alignment.query
score = initial_alignment.score

coordinates = initial_alignment.coordinates
subject_start = coordinates[0,0]
subject_end = coordinates[0,-1]
query_start = coordinates[1,0]
query_end = coordinates[1,-1]

#raw sequence lengths
raw_subject_seq_len = len(subject_seq)
raw_query_seq_len = len(query_seq)

#identify mismatch & gap positions (misalignments)
##match_line_count_chars = ['|', '-', '.']

strip_final_matched_lines = re.sub(r'[^\|\.\-]', '', final_matched_lines)

total_bases = len(strip_final_matched_lines)
base_matches = final_matched_lines.count('|')
misaligned_positions = []

#calc percentage of alignment between subject and query (including gaps)
if total_bases > 0:
    align_percent = round((base_matches / total_bases) * 100, 2)
else:
    align_percent = 'N/A'

#total # of misalignments
total_misalignments = total_bases - base_matches

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
    formatted_alignment=formatted_alignment,
    score=score,
    subject_start=subject_start,
    subject_end=subject_end,
    query_start=query_start,
    query_end=query_end,
    raw_subject_seq_len=raw_subject_seq_len,
    raw_query_seq_len=raw_query_seq_len,
    align_percent=align_percent,
    misaligned_positions=misaligned_positions,
    total_misalignments=total_misalignments,
    final_matched_lines=final_matched_lines,
    pass_status=pass_status,
    highlighted_subject=highlighted_subject,
    highlighted_query=highlighted_query,
    highlighted_match=highlighted_match,
    checkbox_lines=checkbox_lines))
