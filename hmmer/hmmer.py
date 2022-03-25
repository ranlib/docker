#!/usr/bin/env python
'''Run hmmer'''
import os
import argparse
import subprocess
import textwrap
import pyhmmer
import pyrodigal
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def prodigal(infile: str) -> str:
    '''Run prodigal in subprocess.run'''
    name = os.path.splitext(infile)[0]
    faafile = name + '.faa'
    cmd = 'prodigal -q -i {i} -a {n} -p meta'.format(i=infile, n=faafile)
    subprocess.check_output(cmd, shell=True)
    return faafile

def hmmer2(infile: str, hmmer_file: str, threads: int) -> []:
    '''Run hmmer'''
    with pyhmmer.easel.SequenceFile(infile,format='fasta') as file:
        alphabet = file.guess_alphabet()
        #print(alphabet)
        s = pyhmmer.easel.TextSequence()
        sequences = []
        while file.readinto(s) is not None:
            #print(s.sequence)
            ds = s.digitize(alphabet)
            sequences.append(ds)

    print(f"# of sequences {n}".format(n=len(sequences)))
    return len(sequences)

def hmmer(infile: str, hmmer_file: str, threads: int) -> []:
    '''Run hmmer'''
    print(infile)
    with pyhmmer.easel.SequenceFile(infile,format='fasta') as file:
        alphabet = file.guess_alphabet()
        #print(alphabet)
        s = pyhmmer.easel.TextSequence()
        sequences = []
        while file.readinto(s) is not None:
            #print(s.sequence)
            ds = s.digitize(alphabet)
            sequences.append(ds)

    print(f"# of sequences {n}".format(n=len(sequences)))
              
    with pyhmmer.plan7.HMMFile(hmmer_file) as hmms:
        all_hits = list(pyhmmer.hmmsearch(hmms, sequences, cpus=threads))
        #pipeline = pyhmmer.plan7.Pipeline(alphabet)
        #all_hits = [ pipeline.search_hmm(hmm, sequences) for hmm in hmms ]
    return all_hits

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='run hmmer', prog='hmmer', formatter_class=lambda prog: argparse.HelpFormatter(prog,max_help_position=50))
    parser.add_argument('--fasta'   , '-f', type=str, dest='fasta', help='fasta input file', required=True)
    parser.add_argument('--hmm'     , '-m', type=str, dest='hmmer', help='hmmer input file, e.g. AntiFam.hmm', required=True)
    parser.add_argument('--output'  , '-o', type=str, dest='output', help='output file', required=True)
    parser.add_argument('--threads' , '-t', type=int, dest='threads', default=1, help="Number of threads to use. (default: %(default)i)")
    parser.add_argument('--prodigal', '-p', help='run prodigal, requiered if fasta file not protein fasta', action = 'store_true')
    parser.add_argument('--verbose' , '-v', help='modify output verbosity', action = 'store_true')
    args = parser.parse_args()
    if args.prodigal:
        #faa = prodigal(args.fasta)
        with open(args.fasta, "r", encoding="ascii") as fasta:
            gene_records = []
            for record in SeqIO.parse(fasta, "fasta"):
                genes = pyrodigal.Pyrodigal(meta=True).find_genes(str(record.seq))
                for i, gene in enumerate(genes):
                    record = SeqRecord(Seq(gene.translate()), "gene_%i" % (i + 1), "", "")
                    if args.verbose:
                        print(record)
                    gene_records.append(record)
            SeqIO.write(gene_records, args.output, "fasta")
            faa = args.output
    else:
        faa = args.fasta
    hits = hmmer2(faa, args.hmmer, args.threads)
    #print(len(hits))
    # for ihit, hit in enumerate(hits):
    #     for domain in hit.domains:
    #         if domain.score > 9 and len(domain.alignment) > 30:
    #             print(hit.name, domain.alignment.hmm_name)
    #             print(domain.alignment.hmm_sequence)
    #             print(domain.alignment.identity_sequence)
    #             print(domain.alignment.target_sequence)
