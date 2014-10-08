# About 

Command line Blast tutorial written by Yemin Lan (<yeminlan@gmail.com >).

# Notes 

## Countbases in BioPython

```python
#!/usr/bin/env python 
from Bio import SeqIO

sequences = {s.id:s.seq for s in SeqIO.parse("myseqs.fna", "fasta")}

for s_name in sequences.keys(): 
  print s_name + ": " + str(len(sequences[s_name])) + " bases"
``` 
