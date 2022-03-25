#!/usr/bin/env python

from Bio import SeqIO
from enum import Enum

class fastx(Enum):
     FASTA = 1
     FASTQ = 2

def count_bases_fastx(filename: str, filetype: str = "fasta") -> int:
    """count number of bases in fasta or fastq file"""
    return sum( [ len(seq_record) for seq_record in SeqIO.parse(open(filename,'r'), filetype) ] )


def count_reads(filename: str) -> int:
    """count reads"""
    count = 0
    if get_sequence_file_type(filename) == "FASTA":
        for _ in iterate_fasta(filename):
            count += 1
    elif get_sequence_file_type(filename) == "FASTQ":
        for _ in iterate_fastq(filename):
            count += 1
    else:
        sys.exit("<E> count_reads: Error: {} is not FASTA/FASTQ format".format(filename))
    return count


def get_sequence_file_type(filename: str) -> str:
    """
    Determines whether a file is FASTA or FASTQ.
    """
    if not os.path.isfile(filename):
        sys.exit("Error: could not find {}".format(filename))
    if get_compression_type(filename) == "gz":
        open_func = gzip.open
    else:  # plain text
        open_func = open
    with open_func(filename, "rt") as seq_file:
        try:
            first_char = seq_file.read(1)
        except UnicodeDecodeError:
            first_char = ""
    if first_char == ">":
        return "FASTA"
    elif first_char == "@":
        return "FASTQ"
    else:
        raise ValueError("<E> get_sequence_file_type: File is neither FASTA or FASTQ")


def get_open_func(filename: str):
    """test which open function to use"""
    if get_compression_type(filename) == "gz":
        return gzip.open
    else:  # plain text
        return open


def iterate_fasta(filename: str) -> []:
    """iterate fasta file"""
    if get_sequence_file_type(filename) != "FASTA":
        log("<E> iterate_fasta: Error: {} is not FASTA format".format(filename))
        raise StopIteration
    with get_open_func(filename)(filename, "rt") as fasta:
        for line in fasta:
            line = line.strip()
            if len(line) == 0:
                continue
            if not line.startswith(">"):
                continue
            name = line[1:].split()[0]
            sequence = next(fasta).strip()
            yield name, sequence

def iterate_fastq(filename: str) -> []:
    """iterate fastq file"""
    if get_sequence_file_type(filename) != "FASTQ":
        log("<E> iterate_fastq: Error: {} is not FASTQ format".format(filename))
        raise StopIteration
    with get_open_func(filename)(filename, "rt") as fastq:
        for line in fastq:
            line = line.strip()
            if len(line) == 0:
                continue
            if not line.startswith("@"):
                continue
            name = line[1:].split()[0]
            sequence = next(fastq).strip()
            _ = next(fastq)
            qualities = next(fastq).strip()
            yield name, sequence, qualities


