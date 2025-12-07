# FIXALIGN

## About

Web front-end tool for pairwise sequence alignment visualization, approval, and database storage.

Source code is available within this repository.

*Utilizes:*
- Python CGI scripts for backend
- MYSQL database for storage of results
- Jinja2 for HTML templates
- JavaScript for interactive misalignment approval

Alignment and scoring performed using Biopython's PairwiseAligner.

Results can be exported for later use.

## Requirements
- CGI support on webserver (i.e. Apache)
- Python3
- MYSQL database
- Storage space for results data
   * CSS for browser styling
   * Biopython; PairwiseAligner

## Usage

1. If you have not already, download MYSQL & configure your MYSQL credentials within export.cgi; then, create a table in your MYSQL database to hold results variables like so:

```sql CREATE TABLE Alignment_Results (
    id INT AUTO_INCREMENT PRIMARY KEY,
    subject_seq TEXT NOT NULL,
    query_seq TEXT NOT NULL,
    pairwise_score FLOAT,
    scoring VARCHAR(12),
    alignment_score FLOAT,
    pass_status enum('Pass','Fail')
);
```

2. Open input form in browser and enter you subject & query sequences, as well as your alignment parameters. 
3. Submit the form for alignment results viewing.
4. Approve misalignments (mismatches or gaps) by clicking on the misalignment highlighted in red. The percent identity (alignment score) and pass/ fail status will also update automatically.
5. Save results to MYSQL database for downstream usage.
6. View confirmation and saved data results.

For more info on pairwise alignment and scoring, see:
- Biopython documentation https://biopython.org/wiki/Documentation
- General sequence alignment info https://www.ncbi.nlm.nih.gov/BLAST/tutorial/Altschul-1.html

