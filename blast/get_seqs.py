from Bio import Entrez
from Bio import SeqIO
from Bio import Seq
from Bio.Alphabet import IUPAC

Entrez.email = "fake@drexel.edu"
handle = Entrez.esearch(db="nucleotide", term="Escherichia coli K-12[ORGN] AND complete genome")
records = Entrez.read(handle)
handle.close()

handle = Entrez.efetch(db="nucleotide", id=records['IdList'][0], rettype="gb", retmode="text")
record = SeqIO.read(handle, "genbank")
handle.close()

locations=[]
output_handle=open("myseqs.fna","w")
for feature in record.features:
    if feature.type=='CDS':
        if 'gene' in feature.qualifiers:
            if str(feature.location) not in locations:
                locations.append(str(feature.location))
                if len(locations)<=5:
                    output_handle.write(">%s %s\n%s\n" % (feature.qualifiers['gene'][0], record.description, str(feature.extract(record.seq))))
                else:
                    break

output_handle.close()

