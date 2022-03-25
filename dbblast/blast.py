#!/usr/bin/env python
#
# blastn sequence data against nt
#
import sys
import argparse
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastxCommandline
import ssl

ssl._create_default_https_context = ssl._create_unverified_context

#
# input
#
parser = argparse.ArgumentParser(description='blast', prog='blast', formatter_class=lambda prog: argparse.HelpFormatter(prog,max_help_position=50))
parser.add_argument('-i', '--input'   , dest='input'   , type=str, help='Fasta file with sequences', required=True)
parser.add_argument('-o', '--output'  , dest='output'  , type=str, help='Blast results (default=%(default)s)', default="results.xml")
parser.add_argument('-s', '--service' , dest='service' , type=str, help='Type of blast program (default=%(default)s)', default="blastn")
parser.add_argument('-e', '--evalue'  , dest='evalue'  , type=float, help='Evalue (default=%(default)f)', default=1e-20)
parser.add_argument('-d', '--database', dest='database', type=str, help='Database (default=%(default)s)', default='nr')
parser.add_argument('-t', '--threads' , dest='threads' , type=int, help='Number of threads (default=%(default)i)', default=1)
if len(sys.argv)==1:
     parser.print_help(sys.stderr)
     sys.exit(1)
args = parser.parse_args()

with open(args.input) as input_file:
     result_handle = NCBIWWW.qblast(args.service, "nt", input_file.read() ) 

# a different way to do the same
#from Bio import SeqIO 
#seq_record = next(SeqIO.parse(open('test.fasta'),'fasta')) 
#result_handle = NCBIWWW.qblast("blastn", "nt", seq_record.seq) 

#
# write blast result to file
#
with open(args.output, 'w') as save_file: 
     blast_results = result_handle.read() 
     save_file.write(blast_results)

#
# parse blast restuls
#
for record in NCBIXML.parse(open(args.output)): 
     if record.alignments: 
        print("\n") 
        print("query: %s" % record.query[:100]) 
        for align in record.alignments: 
           for hsp in align.hsps: 
              if hsp.expect < args.evalue: 
                 print("match: %s " % align.title[:100])

#
# command line blast
#
blastx_cline = NcbiblastxCommandline(query=args.input, db=args.database, evalue=args.evalue, outfmt=5, out=args.output)
print(blastx_cline)
stdout, stderr = blastx_cline()
