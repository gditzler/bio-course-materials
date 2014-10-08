#!/bin/bash

# module load python/2.7.6
# module load ncbi-blast/gcc/64/2.2.29

# 01.acquire 5 protein-coding genes from a Escherichia coli complete genome and save in myseqs.fna (BioPython)
python get_seqs.py

# 02.check the length of each seq in myseqs.fna (Python)
python countbases.py myseqs.fna

# 03a.blast over the internet using blast+ standalone tools (BinaryExecutable)
blastn -h
blastn -remote -db nr -query myseqs.fna -out myseqs.fna.blastout -evalue 0.00001 -outfmt 6 -max_target_seqs 5
# -perc_identity 97
# -best_hit_overhang 0.25 -best_hit_score_edge 0.25

# 03b.translate nucleotide to protein, blast over the internet using blast+ standalone tools (Python, BinaryExecutable)
# blastp/blastx/tblastn/tblastx
python translate.py myseqs.fna > myseqs.faa
python countbases.py myseqs.faa
blastp -remote -db swissprot -query myseqs.faa -out myseqs.faa.blastout -evalue 0.00001 -outfmt 6 -max_target_seqs 5

# 03c.blast locally using blast+ standalone tools (BinaryExecutable)
blastn -db /path_to_local_db/pdb -query myseqs.fna -out myseqs.fna.blastout -evalue 0.00001 -outfmt 6 -max_target_seqs 5
# -num_threads 1

# 03d.translate nucleotide to protein, blast over the internet using blast+ standalone tools integrated by BioPython (BioPython)
python biopython_blast.py

# 04.blast output formats
# XML vs. tabular

# 05a.acquire only hits mapped to GenBank from blast tabular outputs, acquire the best hit for each query, count the number of queries mapped to each GenBank entry, acquire the GenBank entry info from the gi number (Bash, Eutils)
# parse blastn output myseqs.fna.blastout
grep "|gb|" myseqs.fna.blastout > myseqs.fna.blastout.gbhits
sort -u -k1,1 myseqs.fna.blastout.gbhits > myseqs.fna.blastout.besthit
awk '{print $2}' myseqs.fna.blastout.besthit | sort | uniq -c > myseqs.fna.blastout.besthit.count
for i in `awk -F "|" '{print $2}' myseqs.fna.blastout.besthit.count`; do curl -s "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=$i&retmode=xml" | grep GBSeq_definition | cut -d '>' -f 2 | cut -d '<' -f 1; done
# parse blastp output myseqs.faa.blastout
sort -u -k1,1 myseqs.faa.blastout > myseqs.faa.blastout.besthit
awk '{print $2}' myseqs.faa.blastout.besthit | awk -F "|" '{print $5}'

# 05b.acquire top 3 hits for each query from blast XML outputs, parse for the query name, hit name, evalue, alignment length and the alignment (BioPython)
python biopython_blastparse_xml.py

# 05c.acquire all hits for each query from blast tabular outputs, parse for the query name and hit names (BioPython)
python biopython_blastparse_tabular.py


