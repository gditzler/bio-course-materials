from Bio import SearchIO

blast_records = SearchIO.parse('myseqs.faa.blastout', 'blast-tab')

# parse for top 3 HSPs per query and show the first 75bp of alignment
for blast_record in blast_records:
    print('query: ' + blast_record.id)
    for hit in blast_record.hits:
        print('hit: '+hit.id)



