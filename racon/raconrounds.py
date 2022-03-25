#!/usr/bin/env python3
"""
racon rounds
run minimap2 alignment + racon polishing a number of specified rounds
"""
import os
import sys
import gzip
import argparse
import pathlib
import subprocess
import multiprocessing
import tempfile
import shutil
import platform
import resource
import pysam
import mappy
import racon_utils

def main(args=None):
    """
    run racon
    """
    args = get_arguments(args)
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_dir = pathlib.Path(tmp_dir)
        if args.bam:
            fastq = pysam.fastq(args.bam)
            reads_file = tmp_dir / "temp.fastq"
            with open(reads_file, "wt") as reads:
                number_of_lines = reads.write(fastq)
                log(f"<I> raconrounds: number of lines = {number_of_lines}")
        elif args.reads:
            reads_file = args.reads

        mapping_option = { "hifi": "asm20", "clr": "map-pb", "ont": "map-ont"}
        preset = mapping_option[args.read_type]
        print(args)
        if args.pacbio:
            preset = "map-bp"
        polish(args.assembly, reads_file, args.threads, args.rounds, tmp_dir, preset, args.output, args.trimming, args.include_unpolished, args.debug, args.final)
    print("<I> raconrounds: Peak mem: %s GB" % round(max_mem_usage(), 2))


def max_mem_usage() -> float:
    """Return max mem usage (GB) of self and child processes"""
    max_mem_self = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    max_mem_child = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    if platform.system() == "Linux":
        max_mem = (max_mem_self + max_mem_child) / float(1e6)
    else:
        max_mem = (max_mem_self + max_mem_child) / float(1e9)
    return max_mem


def check_positive(value: int) -> int:
    """check integer is positive"""
    ivalue = int(value)
    if ivalue <= 0:
        raise argparse.ArgumentTypeError("%s is an invalid positive int value" % value)
    return ivalue


def get_default_thread_count() -> int:
    """get default thread count"""
    return min(multiprocessing.cpu_count(), 16)


def get_arguments(args):
    """
    get arguments
    """
    parser = argparse.ArgumentParser(description="raconrounds", prog="raconrounds", formatter_class=lambda prog: argparse.HelpFormatter(prog, max_help_position=80))

    required = parser.add_argument_group("required arguments")
    required.add_argument("-a", "--assembly", dest="assembly", type=str, help="Assembly to be polished (FASTA format)", required=True)
    required.add_argument("-o", "--output", dest="output", type=str, help="Output directory", required=True)

    mutually_exclusive = parser.add_mutually_exclusive_group(required=True)
    mutually_exclusive.add_argument("-r", "--reads", dest="reads", type=str, help="Long reads for polishing (FASTA or FASTQ format)")
    mutually_exclusive.add_argument("-b", "--bam", dest="bam", type=str, help="Long reads for polishing (BAM format)")

    parser.add_argument("-t", "--threads", dest="threads", type=check_positive, help="Number of threads to use for alignment and polishing (default: %(default)i)", default=get_default_thread_count())
    parser.add_argument("-n", "--rounds", dest="rounds", type=check_positive, help="Number of full Racon polishing rounds (default: %(default)i)", default=1)
    parser.add_argument("-s", "--read_type", dest="read_type", type=str, choices=["clr","hifi","ont"], help="Sequencing technology of the reads used (default: %(default)s)", default="hifi")
    parser.add_argument("-p", "--pacbio", action="store_true", help="Use this flag for PacBio reads to make use of the map-pb Minimap2 preset (default=assume Nanopore reads and use the map-ont preset)")
    parser.add_argument("-c", "--trimming", action="store_true", help="turn on racon trimming (default: %(default)s)", default=False)
    parser.add_argument("-i", "--include_unpolished", action="store_true", help="Contigs that could not be polished will be copied unchanged to output (default: %(default)s)", default=True)
    parser.add_argument("-d", "--debug", action="store_true", help="write out intermediary files (default: %(default)s)", default=False)
    parser.add_argument("-f", "--final", action="store_true", help="write out all final data (default: %(default)s)", default=False)

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()
    return args


def polish(assembly_filename: str, read_filename: str, threads: int, rounds: int, tmp_dir: str, preset: str, output: str, trimming: bool = False, include_unpolished: bool = True, debug: bool = False, final: bool = False):
    """
    run racon multiple times
    """
    try:
        os.makedirs(output, exist_ok=True)
    except OSError:
        print("<E> polish: Error creating output directory " + output)
        sys.exit(1)
    i = 0
    polished_filename = assembly_filename
    while i < rounds:
        i = i + 1
        round_name = f"round_{i}"
        unpolished_filename = tmp_dir / (round_name + ".fasta")
        shutil.copyfile(polished_filename, unpolished_filename)
        polished_filename, racon_log, alignments, minimap2_log = run_racon(round_name, read_filename, unpolished_filename, tmp_dir, threads, preset, trimming, include_unpolished)
        if debug:
            for this_file in [polished_filename, racon_log, alignments, minimap2_log]:
                try:
                    shutil.copy(this_file, pathlib.Path(output))
                except shutil.SameFileError:
                    print("<E> polish: Source and destination represents the same file.")
                except PermissionError:
                    print("<E> polish: Permission denied.")

    if not debug and not final:
        for this_file in [polished_filename, racon_log]:
            try:
                shutil.copy(this_file, pathlib.Path(output))
            except shutil.SameFileError:
                print("<E> polish: Source and destination represents the same file.")
            except PermissionError:
                print("<E> polish: Permission denied.")
    if final:
        for this_file in [polished_filename, racon_log, alignments, minimap2_log]:
            try:
                shutil.copy(this_file, pathlib.Path(output))
            except shutil.SameFileError:
                print("<E> polish: Source and destination represents the same file.")
            except PermissionError:
                print("<E> polish: Permission denied.")


def run_racon(name: str, read_filename: str, unpolished_filename: str, tmp_dir: str, threads: int = 1, preset: str = "asm20", trimming: bool = False, include_unpolished: bool = True) -> []:
    """
    run racon
    """
    #read_count = count_reads(read_filename)
    #if read_count <= 1:
    #    log(f"<W> run_racon: Skipping Racon for {name}, number of reads == {read_count} <= 1")
    #    return unpolished_filename

    log(f"<I> run_racon: Running Racon on {name}:")
    #log(f"  reads:      {read_filename} ({read_count:,} reads)")

    unpolished_base_count = count_fasta_bases(unpolished_filename)
    log(f"  input:      {unpolished_filename} ({unpolished_base_count:,} bp)")

    unpolished_base_count = racon_utils.count_bases_fastx(unpolished_filename, "fasta")
    log(f"  input:      {unpolished_filename} ({unpolished_base_count:,} bp)")

    # Align with minimap2
    #preset = "map-pb" if pacbio else "map-ont"
    command = ["minimap2", "-t", str(threads), "-x", preset, "-a", unpolished_filename, read_filename]
    alignments = tmp_dir / (name + ".sam")
    minimap2_log = tmp_dir / (name + "_minimap2.log")
    log("  " + " ".join([str(i) for i in command]))
    with open(alignments, "wt") as stdout, open(minimap2_log, "w") as stderr:
        subprocess.call(command, stdout=stdout, stderr=stderr)

    # Polish with Racon
    polished_filename = tmp_dir / (name + "_polished.fasta")
    command = ["racon", "-t", str(threads)]
    if include_unpolished:
        command.append("--include-unpolished")
    if not trimming:
        command.append("--no-trimming")
    command.extend([read_filename, str(alignments), unpolished_filename])
    racon_log = tmp_dir / (name + "_racon.log")
    log("  " + " ".join([str(i) for i in command]))
    with open(polished_filename, "wt") as stdout, open(racon_log, "w") as stderr:
        subprocess.call(command, stdout=stdout, stderr=stderr)
    polished_base_count = count_fasta_bases(polished_filename)
    log(f"  output:     {polished_filename} ({polished_base_count:,} bp)")
    log()
    return polished_filename, racon_log, alignments, minimap2_log


def log(message: str = "", end: str = "\n"):
    """print message"""
    print(message, file=sys.stderr, flush=True, end=end)


def get_compression_type(filename: str) -> str:
    """
    Attempts to guess the compression (if any) on a file using the first few bytes.
    http://stackoverflow.com/questions/13044562
    """
    magic_dict = {"gz": (b"\x1f", b"\x8b", b"\x08"), "bz2": (b"\x42", b"\x5a", b"\x68"), "zip": (b"\x50", b"\x4b", b"\x03", b"\x04")}
    max_len = max(len(x) for x in magic_dict)

    with open(str(filename), "rb") as unknown_file:
        file_start = unknown_file.read(max_len)

    compression_type = "plain"
    for file_type, magic_bytes in magic_dict.items():
        if file_start.startswith(magic_bytes):
            compression_type = file_type
    if compression_type == "bz2":
        sys.exit("Error: cannot use bzip2 format - use gzip instead")
    if compression_type == "zip":
        sys.exit("Error: cannot use zip format - use gzip instead")
    return compression_type


def load_fasta(fasta_filename: str) -> []:
    """load fasta file"""
    if get_compression_type(fasta_filename) == "gz":
        open_func = gzip.open
    else:  # plain text
        open_func = open
    fasta_seqs = []
    with open_func(fasta_filename, "rt") as fasta_file:
        name = ""
        sequence = []
        for line in fasta_file:
            line = line.strip()
            if not line:
                continue
            if line[0] == ">":  # Header line = start of new contig
                if name:
                    fasta_seqs.append((name.split()[0], "".join(sequence)))
                    sequence = []
                name = line[1:]
            else:
                sequence.append(line)
        if name:
            fasta_seqs.append((name.split()[0], "".join(sequence)))
    return fasta_seqs


def count_fasta_bases(fasta_filename: str) -> int:
    """count number of bases in fasta file"""
    return sum(len(seq) for _, seq in load_fasta(fasta_filename))


if __name__ == "__main__":
    main()
