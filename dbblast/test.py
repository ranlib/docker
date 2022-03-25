
#
# blastn sequence data against nt
#
from Bio.Blast import NCBIWWW
sequence_data = open("test.fasta").read() 
result_handle = NCBIWWW.qblast("blastn", "nt", sequence_data) 

# a different way to do the same
#from Bio import SeqIO 
#seq_record = next(SeqIO.parse(open('test.fasta'),'fasta')) 
#result_handle = NCBIWWW.qblast("blastn", "nt", seq_record.seq) 

#
# write blast result to file
#
with open('results.xml', 'w') as save_file: 
     blast_results = result_handle.read() 
     save_file.write(blast_results)


#
# parse blast restuls
#
from Bio.Blast import NCBIXML
E_VALUE_THRESH = 1e-20 
for record in NCBIXML.parse(open("results.xml")): 
     if record.alignments: 
        print("\n") 
        print("query: %s" % record.query[:100]) 
        for align in record.alignments: 
           for hsp in align.hsps: 
              if hsp.expect < E_VALUE_THRESH: 
                 print("match: %s " % align.title[:100])

