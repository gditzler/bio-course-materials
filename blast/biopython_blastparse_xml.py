from Bio.Blast import NCBIXML

# read blast output
result_handle = open("biopython_blast.xml")
blast_records = NCBIXML.parse(result_handle)
# result_handle.close()

# parse for top 3 HSPs per query and show the first 75bp of alignment
for blast_record in blast_records:
    print('***********************************************')
    print('query: ' + str(blast_record.query))
    for alignment in blast_record.alignments[0:3]:
        for hsp in alignment.hsps:
            print('-----------')
            print('hit: ' + str(alignment.title.split()[0]))
            print('evalue: ' + str(hsp.expect))
            print('length: ' + str(alignment.length))
            if len(hsp.query)>75:
                print('query ' + hsp.query[0:75] + '...')
                print('hit   ' + hsp.sbjct[0:75] + '...')
            else:
                print('query ' + hsp.query)
                print('hit   ' + hsp.sbjct)
            break


