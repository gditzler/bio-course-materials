from Bio.Blast import NCBIWWW
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

# read nucleotides
handle = open("myseqs.fna", "rU")
records = list(SeqIO.parse(handle, "fasta"))
handle.close()

# translate nucleotides to proteins
for i in range(0,len(records)):
    records[i].seq = records[i].seq.translate()

# define output file
save_file = open("biopython_blast.xml", "w")

# blast over the internet
for record in records:
    result_handle = NCBIWWW.qblast("blastp", "swissprot", record, expect=0.00001, format_type='XML', hitlist_size=5)
    save_file.write(result_handle.read())
    result_handle.close()

save_file.close()

