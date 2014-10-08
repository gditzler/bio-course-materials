from Bio import Entrez
from Bio import SeqIO
from Bio import Seq
from Bio.Alphabet import IUPAC


search_terms = [
    "Escherichia coli Nissle 1917, complete genome",
    "Escherichia coli LY180, complete genome",
    "Escherichia coli str. K-12 substr. MC4100 complete genome"
    ] 

for n, term in enumerate(search_terms):
  Entrez.email = "fake@drexel.edu"
  handle = Entrez.esearch(db="nucleotide", term=term)
  records = Entrez.read(handle)
  handle.close()
  
  handle = Entrez.efetch(db="nucleotide", id=records['IdList'][0], rettype="gb", retmode="text")
  record = SeqIO.read(handle, "genbank")
  handle.close()
  
  locations=[]
  output_handle=open("myseqs"+str(n+1)+".fna","w")
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

