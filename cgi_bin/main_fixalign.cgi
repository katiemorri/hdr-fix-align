#!/usr/bin/env python

# FixAlign - Analyze sequence aligments, and easily approve mismatches

import cgi
import jinja2
from Bio import SeqIO
from Bio.pairwise2 import format_alignment

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

#perform alignment via pairwise2, biopython
#gap penalties are preset in this script to -4 and -1 (open, and then penalty)
alignments = pairwise2.align.globalms(subject_seq,
        query_seq,
        match_score,
        mismatch_score,
        -4, -1,
        one_alignment_only=True)

#assign alignment variables (formatted to pairwise output)
initial_alignment = alignments[0]
aligned_subject = initial_alignment.seqA
aligned_query = initial_alignment.seqB
score = initial_alignment.score
start = initial_alignment.start
end = initial_alignment.end

#identify mismatch & gap positions (misalignments)
base_matches = 0
total_bases = len(aligned_subject)
misaligned_positions = []
#get match lines for visualization output between sub and query seq in results form
matched_lines = []

for i in range(total_bases):
    subject_base = aligned_subject[i]
    query_base = aligned_query[i]
    
    #check if not aligned, then save position
    if subject_base != query_base:      #mis-match
        misaligned_positions.append(i)
        matched_lines.append(' ')
    else:
        base_matches += 1               #save for calc. percentage
        matched_lines.append('|')       #match

#determine percentage of alignment between subject and query (including gaps)
align_percent = round((base_matches / total_bases) * 100, 2)

#total # of misalignments
total_misalignments = len(misaligned_positions)

#create line of matches, final
final_matched_lines = ''.join(matched_lines)

#create function to block aligned sub, match lines, and aligned seq
def blocked_lines(seq, block=60):
    return [seq[i:i+block] for i in range(0, len(seq), block)]
blocked_subject = blocked_lines(aligned_subject)
blocked_matched_lines = blocked_lines(final_matched_lines)
blocked_query = blocked_lines(aligned_query)

#determine pass/fail based on % match thresholds
pass_status = "Fail"
if clone == "forward":
    if arm_type == "5Arm":
        if align_percent >= 95:
            pass_status = "Pass"
    elif arm_type == "3Arm":
        if align_percent >= 98:
            pass_status = "Pass"
elif clone == "reverse":
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
    start=start,
    end=end,
    align_percent=align_percent,
    misaligned_positions=misaligned_positions,
    total_misalignments=total_misalignments,
    final_matched_lines=final_matched_lines,
    pass_status=pass_status,
    blocked_subject=blocked_subject,
    blocked_matched_lines=blocked_matched_lines,
    blocked_query=blocked_query))
