#!/usr/bin/env python
from Bio import SeqIO

sequences = SeqIO.parse(open("assembly.faa"),'fasta')
for record in sequences:
    sequence = str(record.seq)
    print(len(sequence))
    
        
