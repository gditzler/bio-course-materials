from Bio import Entrez
from Bio import SeqIO
from Bio import Seq
from Bio.Alphabet import IUPAC

genomes = ["Escherichia coli str. K-12 substr. MC4100 complete genome","Escherichia coli Nissle 1917, complete genome","Escherichia coli LY180, complete genome"]
genomes_short = ["K12","Nissle","LY180"]

for n,genome in enumerate(genomes):
    Entrez.email = "fake@drexel.edu"
    handle = Entrez.esearch(db="nucleotide", term=genome)
    records = Entrez.read(handle)
    handle.close()

    handle = Entrez.efetch(db="nucleotide", id=records['IdList'][0], rettype="gb", retmode="text")
    record = SeqIO.read(handle, "genbank")
    handle.close()

    mygenes = ["thrA","mogA","dnaK","nhaA","ksgA"]
    output_handle=open("seq"+str(n+1)+".fna","w")
    for feature in record.features:
        if feature.type=='CDS':
            if 'gene' in feature.qualifiers:
                if feature.qualifiers['gene'][0] in mygenes:
                    output_handle.write(">%s_%s\n%s\n" % (feature.qualifiers['gene'][0], genomes_short[n], str(feature.extract(record.seq))))

    output_handle.close()

