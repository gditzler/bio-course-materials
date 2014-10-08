#!/bin/bash

mkdir blast_tutorial
cd blast_tutorial

# 01.Local blast pre-installed
module list
module load ncbi-blast/gcc/64/2.2.29
blastn -h
# blastp/blastx/tblastn/tblastx
# makeblastdb -h

# 02.Database to blast against

# 02a.Download existing database from ftp://ftp.ncbi.nlm.nih.gov/blast/db/
mkdir myblastdb
cd myblastdb
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/pataa.tar.gz
tar zxvf pataa.tar.gz
rm pataa.tar.gz
cd ..

# 02b.Build customized database from fasta file

# 02b_1.Recap: acquire some example fasta files using Biopython, learn how to count the base pairs and how to translate DNA to protein using python
cp /mnt/HA/groups/nsftuesGrp/blast/get_seq.py .
python get_seq.py

cp /mnt/HA/groups/nsftuesGrp/blast/countbases.py .
python countbases.py seq1.fna

cp /mnt/HA/groups/nsftuesGrp/blast/translate.py .
python translate.py seq1.fna > seq1.faa
python countbases.py seq1.faa

# 02b_2.Prepare what to blast (query.fna) and what to blast against (reference.fna)
cat seq1.fna seq2.fna > reference.fna
cp seq3.fna query.fna
python translate.py reference.fna > reference.faa
python translate.py query.fna > query.faa

# 02b_3.Build customized database from fasta files
makeblastdb -h
makeblastdb -in reference.fna -dbtype nucl -out myblastdb/Ecoli
makeblastdb -in reference.faa -dbtype prot -out myblastdb/Ecoli
# "-parse_seqids" if the sequence headers were in NCBI style and therefore can be parsed. For example:
# ">gi|568336023|gb|CM000663.2| Homo sapiens chromosome 1, GRCh38 reference primary assembly."

# 03.Blast
blastn -db ./myblastdb/Ecoli -evalue 0.00001 -outfmt 6 -max_target_seqs 3 -query query.fna -out query.fna.blast1
blastp -db ./myblastdb/Ecoli -evalue 0.00001 -outfmt 6 -max_target_seqs 3 -query query.faa -out query.faa.blast2
blastp -db ./myblastdb/pataa -evalue 0.00001 -outfmt 6 -max_target_seqs 3 -num_threads 8 -query query.faa -out query.faa.blast3
blastn -remote -db nr -evalue 0.00001 -outfmt 6 -max_target_seqs 3 -query query.fna -out query.fna.blast4
# -perc_identity 97
# -best_hit_overhang 0.25 -best_hit_score_edge 0.25
# NOTE: use "-num_threads $NSLOTS" when submitting a proteus job script

# 04.Parse the results, using query.faa.blast2 as an example

# parse for hits matched to E.coli K-12 strain
grep "K12" query.faa.blast2

# parse for best hits for each query
sort -u -k1,1 query.faa.blast2 > query.faa.blast2.besthit

# count the hits
awk '{print $2}' query.faa.blast2.besthit | awk -F "_" '{print $2}' | sort | uniq -c

# parse for hits with >1000 Score
awk '$12 > 1000 { print $0 }' query.faa.blast2 

# 05.A complete example
cat seq*.fna > all.fna
python translate.py all.fna > all.faa
blastp -remote -db swissprot -evalue 0.00001 -outfmt 6 -max_target_seqs 3 -query all.faa -out all.faa.blast5
sort -u -k1,1 all.faa.blast5 > all.faa.blast5.besthit
awk '{print $2}' all.faa.blast5.besthit | awk -F "|" '{print $5}' | sort | uniq -c

# 06.Getting the lowest common ancestor
grep "nhaA_K12" all.faa.blast5

grep "nhaA_K12" all.faa.blast5 | awk -F "|" '{print $2}' > GI.txt
for i in `cat GI.txt`; do curl -s "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=$i&retmode=xml" | grep GBSeq_taxonomy | cut -d '>' -f 2 | cut -d '<' -f 1; done

grep "thrA_Nissle" all.faa.blast5 | awk -F "|" '{print $2}' > GI.txt
for i in `cat GI.txt`; do curl -s "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=$i&retmode=xml" | grep GBSeq_taxonomy | cut -d '>' -f 2 | cut -d '<' -f 1; done

# NCBI Taxa ID
for i in `cat GI.txt`; do curl -s "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=$i&retmode=xml" | grep taxon: | cut -d '>' -f 2 | cut -d '<' -f 1; done

# "Megan5" lowest common ancestor algorithm
# "Blast2lca" on github (written in Go language)
# "Classify BLAST" on Galaxy Tool Shed
# UniProt/Swissprot website http://www.uniprot.org/uploadlists/


