#!/usr/bin/env python
"""
create coverage plots for a given bam file
bam file will be sorted and indexed if not yet sorted
"""
from __future__ import annotations
import sys
import os
import shutil
import argparse
import re
import math
import csv
import platform
import resource
import multiprocessing
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import yaml
import pysam

def check_bam(bam_file: str, nthreads: int = 1, memory: str = "10G", make_new_index: bool = False, verbose: bool = False) -> str:
    """
    check if bam file sorted and index, sort and index bam file if needed
    """
    # check if sorted
    with pysam.AlignmentFile(bam_file, "rb") as aln:
        if verbose:
            log(aln.header)
        if "HD" not in aln.header or ("HD" in aln.header and aln.header["HD"]["SO"] != "coordinate"):
            log("<I> check_bam: sorting bam file " + bam_file)
            sorted_bam = os.path.splitext(bam_file)[0] + ".sorted.bam"
            unsorted_bam = os.path.splitext(bam_file)[0] + ".unsorted.bam"
            pysam.sort("-@", str(nthreads), "-m", memory, "-o", sorted_bam, bam_file)
            shutil.copy(bam_file, unsorted_bam)
            os.rename(sorted_bam, bam_file)

        # if "{}.bai".format(bam_file) not in os.listdir(os.path.dirname(bam_file)) or make_new_index:
        #    print("<I> check_bam: indexing bam file " + bam_file)
        #    pysam.index(bam_file)
    # check if indexed
    # read new alignment object
    if make_new_index:
        log("<I> check_bam: force indexing bam file " + bam_file)
        pysam.index(bam_file)
    else:
        with pysam.AlignmentFile(bam_file, "rb") as alnm:
            try:
                alnm.check_index()
            except ValueError:
                log("<I> check_bam: indexing bam file " + bam_file)
                pysam.index(bam_file)

    return bam_file


def sam_to_bam(samfile: str, nthreads: int = 1, memory: str = "10G") -> str:
    """Convert sam to bam, sort and index bam file"""
    with pysam.AlignmentFile(samfile, "r") as infile:
        unsorted_bam = os.path.splitext(samfile)[0] + ".unsorted.bam"
        with pysam.AlignmentFile(unsorted_bam, "wb", template=infile) as bamfile:
            for line in infile:
                bamfile.write(line)
        sorted_bam = os.path.splitext(samfile)[0] + ".bam"
        pysam.sort("-@", str(nthreads), "-m", memory, "-o", sorted_bam, unsorted_bam)
        pysam.index(sorted_bam)
    return sorted_bam


def get_chromosomes(bam_file: str, keep_contigs: str = None) -> dict[str, int]:
    """Create dictionary of chromosome lengths"""
    with pysam.AlignmentFile(bam_file, "rb") as aln:
        idxstats = aln.get_index_statistics()

        if keep_contigs is None:
            keep_contigs = "."
            pattern = re.compile(keep_contigs)  # only keep contigs that match a certain pattern

        # only keep contigs that match a certain pattern
        # and have mapped reads
        keep_contigs_list : list = []
        for i in idxstats:
            if i.mapped > 0 and pattern.match(i.contig):
                keep_contigs_list.append(i.contig)

        conlen : dict[str,int] = {x: aln.get_reference_length(x) for x in keep_contigs_list}
    return conlen


def get_stats(coverage: dict[str, int], lower: int, upper: int) -> dict[str, str | int | float ]:
    """Get stats, % below and above threshold, for a contig"""
    stats : dict[str, str | int | float] = {"name": "", "maxi": 0, "mini": 0, "mean": 0, "std": 0, "len": 0, "percent_below_threshold": 0, "percent_above_threshold": 0}
    if len(coverage) > 0:
        contig = list(coverage.keys())[0]
        stats["name"] = contig
        stats["maxi"] = int(coverage[contig].max())
        stats["mini"] = int(coverage[contig].min())
        stats["mean"] = round(float(coverage[contig].mean()), 2)
        stats["std"] = round(float(coverage[contig].std()), 2)
        stats["len"] = len(coverage[contig])
        # collect data below and above threshold
        lower_coverage = (coverage[contig])[coverage[contig] < lower]
        upper_coverage = (coverage[contig])[coverage[contig] > upper]
        if len(coverage[contig]) > 0:
            stats["percent_below_threshold"] = round(len(lower_coverage) * 100.0 / len(coverage[contig]), 2)
            stats["percent_above_threshold"] = round(len(upper_coverage) * 100.0 / len(coverage[contig]), 2)
        else:
            stats["percent_below_threshold"] = -1.0
            stats["percent_above_threshold"] = -1.0
    return stats


def get_coverage(bam_file: str, outdir: str, base_quality_cutoff: int, width: int, gwidth: int, lower: int, upper: int, verbose: bool) -> dict[str,float]:
    """Make coverage plots"""
    contigs_data = {}
    coverage_total = {"total": np.array([])}
    with pysam.AlignmentFile(bam_file, "rb") as aln:
        contigs = get_chromosomes(bam_file)
        for contig in contigs:
            coverage = {contig: aln.count_coverage(contig, None, None, None, base_quality_cutoff)}  # Note: coverage = array of 4 arrays, one per A, C, G, T
            coverage_sum = {contig: np.sum(coverage[contig], axis=0)}  # get coverage for sum of A, C, G, T
            coverage_total["total"] = np.append(coverage_total["total"], coverage_sum[contig])
            contigs_data[contig] = get_stats(coverage_sum, lower, upper)

            if verbose:
                # log("<I> bam_coverage: " + contig + " " + str(len(coverage_sum[contig])))
                log("<I> bam_coverage: {name} {len} {mini} {maxi} {mean} {std}".format(name=contigs_data[contig]["name"], len=contigs_data[contig]["len"], mini=contigs_data[contig]["mini"], maxi=contigs_data[contig]["maxi"], mean=contigs_data[contig]["mean"], std=contigs_data[contig]["std"]))

            fig = plt.figure()
            nbins = max(1, math.ceil((coverage_sum[contig].max() - coverage_sum[contig].min()) / width))  # set bin width to width
            plt.hist(coverage_sum[contig], alpha=0.9, color="orange", bins=nbins)
            plt.title(contig + " Coverage")
            plt.xlabel("Number of reads")
            plt.ylabel("Frequency")
            plt.savefig(os.path.join(outdir, contig + "_coverage.png"), bbox_inches="tight", dpi=600)
            plt.close(fig)

            fig = plt.figure()
            plt.title(contig + " Coverage")
            sns_plot = sns.lineplot(x=range(len(coverage_sum[contig])), y=coverage_sum[contig])
            sns_plot.set(xlabel="Genome Position (bp)", ylabel="Number of reads")
            plt.savefig(os.path.join(outdir, contig + "_coverage_vs_loc.png"), bbox_inches="tight", dpi=600)
            plt.close(fig)

            fig = plt.figure()
            data_frame = pd.DataFrame(
                {
                    "X": range(len(coverage_sum[contig])),
                    "Y": coverage_sum[contig],
                }
            )
            mbins = max(1, math.ceil(len(coverage_sum[contig]) / gwidth))  # set bin width to width
            data_frame["Xbins"] = np.digitize(data_frame.X, np.linspace(0, len(coverage_sum[contig]), mbins))
            data_frame["Ymean"] = data_frame.groupby("Xbins").Y.transform("mean")
            axis = data_frame.plot(kind="line", x="Xbins", y="Ymean", title=contig + " Coverage", legend=False)
            axis.set_xlabel("Genome Position (bp)")
            axis.set_ylabel("Number of reads")
            yabs_max = abs(max(axis.get_ylim(), key=abs))
            axis.set_ylim(ymin=-1, ymax=yabs_max)
            plt.savefig(os.path.join(outdir, contig + "_coverage_vs_location_binned.png"), bbox_inches="tight", dpi=600)
            plt.close(fig)
            plt.close("all")

    fig = plt.figure()
    nbins = max(1, math.ceil((coverage_total["total"].max() - coverage_total["total"].min()) / width))  # set bin width to width
    plt.hist(coverage_total["total"], alpha=0.9, color="cornflowerblue", bins=nbins)
    plt.title("Coverage")
    plt.xlabel("Number of reads")
    plt.ylabel("Frequency")
    plt.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    plt.savefig(os.path.join(outdir, "coverage.png"), bbox_inches="tight", dpi=600)
    plt.close(fig)
    plt.close("all")
    # store information
    data = {}
    data["total"] = get_stats(coverage_total, lower, upper)
    data["contigs"] = contigs_data
    return data


def get_coverage_total(bam_file: str, outdir: str, base_quality_cutoff: int, width: int, gwidth: int, lower: int, upper: int) -> dict[str, str | float ]:
    """Make coverage plots"""
    coverage_total = {"total": np.array([])}
    contigs = get_chromosomes(bam_file)
    for contig in contigs:
        with pysam.AlignmentFile(bam_file, "rb") as aln:
            coverage = {contig: aln.count_coverage(contig, None, None, None, base_quality_cutoff)}  # Note: coverage = array of 4 arrays, one per A, C, G, T
            coverage_sum = {contig: np.sum(coverage[contig], axis=0)}  # get coverage for sum of A, C, G, T
            coverage_total["total"] = np.append(coverage_total["total"], coverage_sum[contig])

    fig = plt.figure()
    nbins = max(1, math.ceil((coverage_total["total"].max() - coverage_total["total"].min()) / width))  # set bin width to width
    plt.hist(coverage_total["total"], alpha=0.9, color="cornflowerblue", bins=nbins)
    plt.title("Coverage")
    plt.xlabel("Number of reads")
    plt.ylabel("Frequency")
    plt.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    plt.savefig(os.path.join(outdir, "coverage.png"), bbox_inches="tight", dpi=600)
    plt.close(fig)

    fig = plt.figure()
    data_frame = pd.DataFrame(
        {
            "X": range(len(coverage_total["total"])),
            "Y": coverage_total["total"],
        }
    )
    mbins = max(1, math.ceil(len(coverage_total["total"]) / gwidth))  # set bin width to width
    data_frame["Xbins"] = np.digitize(data_frame.X, np.linspace(0, len(coverage_total["total"]), mbins))
    data_frame["Ymean"] = data_frame.groupby("Xbins").Y.transform("mean")
    axis = data_frame.plot(kind="line", x="Xbins", y="Ymean", title="Coverage", legend=False)
    axis.set_xlabel("Genome Position (bp)")
    axis.set_ylabel("Number of reads")
    yabs_max = abs(max(axis.get_ylim(), key=abs))
    axis.set_ylim(ymin=-1, ymax=yabs_max)
    plt.savefig(os.path.join(outdir, "coverage_vs_location_binned.png"), bbox_inches="tight", dpi=600)
    plt.close(fig)
    plt.close("all")
    return get_stats(coverage_total, lower, upper)


def process_contig(contig: str, bam_file: str, outdir: str, base_quality_cutoff: int, width: int, gwidth: int, lower: int, upper: int) -> dict[str,float]:
    """process contig"""
    contigs_data = {}
    with pysam.AlignmentFile(bam_file, "rb") as aln:
        coverage = {contig: aln.count_coverage(contig, None, None, None, base_quality_cutoff)}  # Note: coverage = array of 4 arrays, one per A, C, G, T
        coverage_sum = {contig: np.sum(coverage[contig], axis=0)}  # get coverage for sum of A, C, G, T
        contigs_data[contig] = get_stats(coverage_sum, lower, upper)

        fig = plt.figure()
        nbins = max(1, math.ceil((coverage_sum[contig].max() - coverage_sum[contig].min()) / width))  # set bin width to width
        plt.hist(coverage_sum[contig], alpha=0.9, color="orange", bins=nbins)
        plt.title(contig + " Coverage")
        plt.xlabel("Number of reads")
        plt.ylabel("Frequency")
        plt.savefig(os.path.join(outdir, contig + "_coverage.png"), bbox_inches="tight", dpi=600)
        plt.close(fig)

        fig = plt.figure()
        plt.title(contig + " Coverage")
        sns_plot = sns.lineplot(x=range(len(coverage_sum[contig])), y=coverage_sum[contig])
        sns_plot.set(xlabel="Genome Position (bp)", ylabel="Number of reads")
        plt.savefig(os.path.join(outdir, contig + "_coverage_vs_loc.png"), bbox_inches="tight", dpi=600)
        plt.close(fig)

        fig = plt.figure()
        data_frame = pd.DataFrame(
            {
                "X": range(len(coverage_sum[contig])),
                "Y": coverage_sum[contig],
            }
        )
        mbins = max(1, math.ceil(len(coverage_sum[contig]) / gwidth))  # set bin width to width
        data_frame["Xbins"] = np.digitize(data_frame.X, np.linspace(0, len(coverage_sum[contig]), mbins))
        data_frame["Ymean"] = data_frame.groupby("Xbins").Y.transform("mean")
        axis = data_frame.plot(kind="line", x="Xbins", y="Ymean", title=contig + " Coverage", legend=False)
        axis.set_xlabel("Genome Position (bp)")
        axis.set_ylabel("Number of reads")
        yabs_max = abs(max(axis.get_ylim(), key=abs))
        axis.set_ylim(ymin=-1, ymax=yabs_max)
        plt.savefig(os.path.join(outdir, contig + "_coverage_vs_location_binned.png"), bbox_inches="tight", dpi=600)
        plt.close(fig)
        plt.close("all")
    return contigs_data


def log(message: str = "", end: str = "\n"):
    """print log message"""
    print(message, file=sys.stderr, flush=True, end=end)


def check_positive(value: int) -> int:
    """check integer is positive"""
    ivalue = int(value)
    if ivalue <= 0:
        raise argparse.ArgumentTypeError("%s is an invalid positive int value" % value)
    return ivalue


def get_default_thread_count() -> int:
    """get default thread count"""
    return min(multiprocessing.cpu_count(), 16)


def max_mem_usage() -> float:
    """Return max mem usage (GB) of self and child processes"""
    max_mem_self = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    max_mem_child = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    if platform.system() == "Linux":
        max_mem = (max_mem_self + max_mem_child) / float(1e6)
    else:
        max_mem = (max_mem_self + max_mem_child) / float(1e9)
    return max_mem


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="make coverage plots", prog=os.path.basename(__file__), formatter_class=lambda prog: argparse.HelpFormatter(prog, max_help_position=80))
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--bam", "-b", type=str, dest="bam_file", help="bam file name")
    group.add_argument("--sam", "-s", type=str, dest="sam_file", help="sam file name")
    parser.add_argument("--output", "-o", type=str, dest="out_dir", help="output directory", required=True)
    parser.add_argument("--upper", "-u", type=check_positive, dest="upper", help="upper threshold for coverage (default: %(default)i)", default=100000000)
    parser.add_argument("--lower", "-l", type=check_positive, dest="lower", help="lower threshold for coverage (default: %(default)i)", default=0)
    parser.add_argument("--width", "-w", type=check_positive, dest="width", help="bin width for read distribution histogram (default: %(default)i)", default=1)
    parser.add_argument("--gwidth", "-g", type=check_positive, dest="gwidth", help="bin width for genome location (default: %(default)i)", default=1)
    parser.add_argument("--threads", "-t", type=check_positive, dest="threads", help="number of threads to use (default: %(default)i)", default=1)
    parser.add_argument("--qual", "-q", type=int, choices=range(0, 100), dest="base_qual", help="minimum base quality to count a base (default: %(default)i)", default=0, metavar="QUAL")
    parser.add_argument("--memory", "-m", type=str, dest="memory", help="set memory for samtools threads (default: %(default)s)", default="10G")
    parser.add_argument("--index", "-i", dest="index", help="force creation of new bam index (default: %(default)s)", action="store_true")
    parser.add_argument("--verbose", "-v", dest="verbose", help="modify output verbosity (default: %(default)s)", action="store_true")
    parser.add_argument("--yaml", "-y", dest="yaml_output", help="enable output of results in yaml format (default: %(default)s)", action="store_true")
    parser.add_argument("--multi_processing", "-p", type=int, choices=range(1, 1000), dest="number_of_cpus", help="enable parallel processing of contigs (default: %(default)i)", default=0, metavar="NCPUS")
    args = parser.parse_args()

    if args.lower > args.upper:
        print("<E> bam_coverage: upper threshold < lower threshold")
        sys.exit(1)

    if args.threads > multiprocessing.cpu_count():
        threads = multiprocessing.cpu_count()
        print("<W> bam_coverage: reset number of threads to multiprocessing CPU count")
    else:
        threads = args.threads

    if not os.path.isdir(args.out_dir):
        print("<I> bam_coverage: Cannot find directory %s, creating as new" % args.out_dir)
        try:
            os.makedirs(args.out_dir, exist_ok=True)
        except OSError:
            print("<E> bam_coverage: Error creating output directory " + args.out_dir)
            sys.exit(1)

    if args.bam_file:  # still need to sort + index
        if os.path.isfile(args.bam_file):
            bam = check_bam(args.bam_file, threads, args.memory, args.index, args.verbose)
        else:
            print("<E> bam_coverage: " + args.bam_file + " does not exist!")
            sys.exit(1)

    if args.sam_file:  # still need sam -> bam + bam sort + index
        if os.path.isfile(args.sam_file):
            bam = sam_to_bam(args.sam_file, threads, args.memory)
        else:
            print("<E> bam_coverage: " + args.sam_file + " does not exist!")
            sys.exit(1)

    if args.number_of_cpus > 0:
        # multiprocessing
        ncontigs = get_chromosomes(bam)
        inputs = [(contig, bam, args.out_dir, args.base_qual, args.width, args.gwidth, args.lower, args.upper) for contig in ncontigs]
        NPROCESSES = min(multiprocessing.cpu_count(), args.number_of_cpus)
        with multiprocessing.Pool(NPROCESSES) as p:
            outputs = p.starmap(process_contig, inputs)
        coverage_stats = {}
        coverage_stats["total"] = get_coverage_total(bam, args.out_dir, args.base_qual, args.width, args.gwidth, args.lower, args.upper)
        coverage_stats["contigs"] = {}
        for contig_dict in outputs:
            coverage_stats["contigs"].update(contig_dict)
    else:
        # loop over contigs
        coverage_stats = get_coverage(bam, args.out_dir, args.base_qual, args.width, args.gwidth, args.lower, args.upper, args.verbose)

    # store information
    if args.yaml_output:
        try:
            with open(os.path.join(args.out_dir, "coverage.yaml"), "wt", encoding="ascii") as yaml_file:
                yaml.dump(coverage_stats, yaml_file, default_flow_style=False, sort_keys=False)
        except IOError:
            print("<E> bam_coverage: I/O error while writing yaml file")

    try:
        with open(os.path.join(args.out_dir, "coverage.tsv"), "w", encoding="ascii") as tsv_file:
            writer = csv.DictWriter(tsv_file, fieldnames=coverage_stats["total"].keys(), delimiter="\t")
            writer.writeheader()
            writer.writerow(coverage_stats["total"])
            for record in coverage_stats["contigs"]:
                writer.writerow(coverage_stats["contigs"][record])
    except IOError:
        print("<E> bam_coverage: I/O error while writing tsv file")

    print("<I> bam_coverage: Peak mem: %s GB" % round(max_mem_usage(), 2))
