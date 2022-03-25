#!/usr/bin/env python
"""
generate plant QC report
"""
import sys
import os
import argparse
from collections import OrderedDict
import locale
import json
import logging
import configparser
from typing import Dict
import nob
import pandas
import pdfkit
import requests
from jira import JIRA
from jira.exceptions import JIRAError

DESCRIPTION = """create table of metrics for pacbio or illumina libraries"""


def get_arguments():
    """get arguments"""
    parser = argparse.ArgumentParser(description=DESCRIPTION, prog=os.path.basename(__file__), formatter_class=lambda prog: argparse.HelpFormatter(prog, max_help_position=80))
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--file", "-f", dest="input_file", type=str, help="text file with library names, one per line")
    group.add_argument("--lib", "-l", dest="library_name", type=str, help="name of library")
    parser.add_argument("--output", "-o", dest="output", type=str, help="output filename (default: %(default)s)", default="plantqc.tsv")
    parser.add_argument("--verbose", "-v", dest="verbose", help="modify output verbosity (default: %(default)s)", action="store_true")
    parser.add_argument("--grouping", "-g", dest="grouping", help="format numbers (default: %(default)s)", action="store_false")
    parser.add_argument("--transpose", "-t", dest="transpose", help="transpose output file (default: %(default)s)", action="store_false")
    parser.add_argument("--excel", "-e", dest="excel", help="add excel file output (default: %(default)s)", action="store_true")
    parser.add_argument("--markdown", "-m", dest="markdown", help="add markdown file output (default: %(default)s)", action="store_false")
    parser.add_argument("--html", "-w", dest="html", help="add html file output (default: %(default)s)", action="store_true")
    parser.add_argument("--pdf", "-p", dest="pdf", help="add pdf file output (default: %(default)s)", action="store_true")
    parser.add_argument("--seq_unit", "-u", dest="seq_unit", help="file output per seq_unit (default: %(default)s)", action="store_true")
    parser.add_argument("--ccs", "-c", dest="ccs", type=int, choices=range(0, 3), help="for pacbio libraries: only ccs (0), subreads (1), both (3) (default: %(default)i)", default=0)
    parser.add_argument("--jira", "-j", dest="jira_ticket", type=str, help="name of jira ticket", default=argparse.SUPPRESS, required=False)
    parser.add_argument("--user", dest="jira_user", type=str, help="jira user name", default=argparse.SUPPRESS, required=False)
    parser.add_argument("--password", dest="jira_password", type=str, help="jira user password", default=argparse.SUPPRESS, required=False)
    parser.add_argument("--log", dest="log", type=str, choices=["debug", "critical", "error", "warning", "info"], help="set log level (default: %(default)s)", default="error", required=False)
    parser.add_argument("--config", dest="config", type=str, help="Configuration file", required=False)
    return parser.parse_args()


def family_of_top_sketch_hit(file_name: str) -> str:
    """get family of top sketch hit"""
    url = file_name.replace("/global/dna/shared/rqc/archive", "https://rqc.jgi.lbl.gov/api/file/file/dna")
    dataframe = pandas.read_csv(url, sep="\t", skiprows=3)
    dataframe = dataframe.dropna()
    d_sorted = dataframe.sort_values(by="Matches", ascending=False)  # sort by number of matches, most matches at the top
    top_tax = d_sorted["taxonomy"][0].split(";")  # take top hit
    result = [record.split(":")[1] for record in top_tax if record.startswith("f:")]
    if len(result) == 0:
        family = "None"
    elif len(result) == 1:
        family = result[0]
    else:
        print("<W> parse_family_from_sketch: more than 1 family record, take only first one.")
        family = result[0]
    return family


def taxonomy_of_top_sketch_hit(file_name: str) -> Dict[str, str]:
    """get family of top sketch hit"""
    url = file_name.replace("/global/dna/shared/rqc/archive", "https://rqc.jgi.lbl.gov/api/file/file/dna")
    dataframe = pandas.read_csv(url, sep="\t", skiprows=3)
    dataframe = dataframe.dropna()
    d_sorted = dataframe.sort_values(by="Matches", ascending=False)  # sort by number of matches, most matches at the top
    return {k.split(":")[0]: k.split(":")[1] for k in d_sorted["taxonomy"][0].split(";")}


def get_pacbio_stats(library_name: str) -> Dict[str, str]:
    """get pacbio statistics"""
    pbstats = "https://rqc.jgi-psf.org/api/rqcws/pbstats/" + library_name
    data = {}
    try:
        response = requests.get(pbstats)
    except requests.exceptions.RequestException as error:
        print("<E> get_pacbio_stats" + str(error))
    else:
        data = response.json()
    return data


def get_pacbio_stats_hash(library_name: str) -> Dict[str, str]:
    """get pacbio statistics"""
    pbstats = "https://rqc.jgi-psf.org/api/rqcws/pbstats/" + library_name
    data = {}
    try:
        response = requests.get(pbstats)
    except requests.exceptions.RequestException as error:
        print("<E> get_pacbio_stats" + str(error))
    else:
        data = response.json()

    # instead of array hash on seq_unit_name
    data_hash = {}
    for seq_unit in data["seq_units"]:
        data_hash[seq_unit["Seq Unit Name"]] = seq_unit

    return data_hash


def get_library_info(library_name: str) -> Dict[str, str]:
    """get library info"""
    data = {}
    library_info = "https://rqc.jgi-psf.org/api/library/info/" + library_name
    logger.debug(library_info)
    try:
        response = requests.get(library_info)
    except requests.exceptions.RequestException as error:
        print("<E> get_library_info: " + str(error))
    else:
        data = response.json()
    return data


def get_report(library_name: str) -> Dict[str, str]:
    """get library report"""
    report = "https://rqc.jgi-psf.org/api/fullreport/report/" + library_name
    logger.debug(report)
    data = {}
    try:
        response = requests.get(report)
    except requests.exceptions.RequestException as error:
        print("<E> get_report: " + str(error))
    else:
        data = response.json()
    return data


def get_meta_data(library_name: str, lims_id: str) -> Dict[str, str]:
    """get meta information"""
    metadata = f"http://ws-access.jgi.doe.gov/pps-pru-library-metadata?pru-lims-id={lims_id}&library-name={library_name}"
    data = {}
    try:
        response = requests.get(metadata)
    except requests.exceptions.RequestException as error:
        print("<E> get_meta_data: " + str(error))
    else:
        data = response.json()
    return data


def pacbio_metrics(lib: str, ccs: int) -> pandas.DataFrame:
    """query library in database"""
    mydict = OrderedDict()

    # get info via rest calls
    data_hash = get_pacbio_stats_hash(lib)
    lib_info = get_library_info(lib)

    # filter lib_info for ccs or subreads
    my_seq_units = []

    if ccs not in [0, 1, 2]:
        print("<W> pacbio_metrics: invalid value for ccs, ccs needs to be in [0, 1, 2]")

    for seq_unit in lib_info["library_info"]:
        seq_unit_name = seq_unit["seq_unit_name"]

        if (ccs == 0) and ("ccs" in seq_unit_name):  # only ccs reads
            my_seq_units.append(seq_unit)

        if (ccs == 1) and ("subreads" in seq_unit_name):  # only subreads
            my_seq_units.append(seq_unit)

        if ccs == 2:
            my_seq_units.append(seq_unit)

    # loop over sequnits
    for seq_unit in my_seq_units:
        seq_unit_name = seq_unit["seq_unit_name"]
        mydict[seq_unit_name] = {}

        mykeys = ["Library Name", "Seq Unit Name", "Sample Family", "Subread Counts", "P0", "P1", "P2", "PostFilter Read Count", "Reads of Insert Length", "File1", "File2"]
        for key in mykeys:
            mydict[seq_unit_name][key] = data_hash[seq_unit_name][key]

        # calculate some derived quantities
        mydict[seq_unit_name]["P0[%]"] = -1
        mydict[seq_unit_name]["P1[%]"] = -1
        mydict[seq_unit_name]["P2[%]"] = -1
        mydict[seq_unit_name]["P0+P1+P2"] = mydict[seq_unit_name]["P0"] + mydict[seq_unit_name]["P1"] + mydict[seq_unit_name]["P2"]
        if mydict[seq_unit_name]["P0+P1+P2"] > 0:
            mydict[seq_unit_name]["P0[%]"] = round(mydict[seq_unit_name]["P0"] * 100 / mydict[seq_unit_name]["P0+P1+P2"], 2)
            mydict[seq_unit_name]["P1[%]"] = round(mydict[seq_unit_name]["P1"] * 100 / mydict[seq_unit_name]["P0+P1+P2"], 2)
            mydict[seq_unit_name]["P2[%]"] = round(mydict[seq_unit_name]["P2"] * 100 / mydict[seq_unit_name]["P0+P1+P2"], 2)

        # check out family
        mydict[seq_unit_name]["Family"] = family_of_top_sketch_hit(mydict[seq_unit_name]["File2"])
        if mydict[seq_unit_name]["Sample Family"] == "None":
            mydict[seq_unit_name]["Correct Family"] = "No Info"
        else:
            mydict[seq_unit_name]["Correct Family"] = "Yes" if mydict[seq_unit_name]["Sample Family"] == mydict[seq_unit_name]["Family"] else "No"

        # add taxonomy & raw-base_count
        fields = ["tax_family", "tax_genus", "tax_species", "gls_physical_run_unit_id", "sdm_raw_base_count"]
        for field in fields:
            mydict[seq_unit_name][field] = seq_unit[field]

        mydict[seq_unit_name]["Raw base count [Gb]"] = round(float(seq_unit["sdm_raw_base_count"]) / 1e9, 1)

        # add logial amounts & genome size
        meta = get_meta_data(lib, mydict[seq_unit_name]["gls_physical_run_unit_id"])
        fields = ["completed-logical-amount", "target-logical-amount", "logical-amount-units", "genome-size-estimated-mb", "sow-item-type"]
        for field in fields:
            if field in meta["libraries-metadata"][0]:
                mydict[seq_unit_name][field] = meta["libraries-metadata"][0][field]

        # check out coverage
        mydict[seq_unit_name]["complete_logical_amount (CLA)"] = round(float(mydict[seq_unit_name]["sdm_raw_base_count"]) / 1e6 / float(mydict[seq_unit_name]["genome-size-estimated-mb"]), 2)

    # create pandas data frame
    dataframe = pandas.DataFrame.from_dict(mydict, orient="index")
    dataframe.drop(columns=["File1", "File2"], inplace=True)
    return dataframe


def illumina_metrics(lib: str, configuration: configparser.ConfigParser) -> pandas.DataFrame:
    """query library in database"""
    mydict = OrderedDict()
    lib_info = get_library_info(lib)
    logger.debug(json.dumps(lib_info, indent=3))
    for seq_unit in lib_info["library_info"]:
        seq_unit_name = seq_unit["seq_unit_name"]
        mydict[seq_unit_name] = {}
        mydict[seq_unit_name]["sdm_raw_base_count"] = seq_unit["sdm_raw_base_count"]
        mydict[seq_unit_name]["Raw base count [Gb]"] = round(float(seq_unit["sdm_raw_base_count"]) / 1e9, 1)
        mydict[seq_unit_name]["gls_physical_run_unit_id"] = seq_unit["gls_physical_run_unit_id"]
        meta = get_meta_data(lib, seq_unit["gls_physical_run_unit_id"])
        logger.debug(json.dumps(meta, indent=3))
        fields = ["completed-logical-amount", "target-logical-amount", "genome-size-estimated-mb", "sow-item-type"]
        for field in fields:
            if field in meta["libraries-metadata"][0]:
                mydict[seq_unit_name][field] = meta["libraries-metadata"][0][field]

    for seq_unit_name in mydict:
        logger.debug(seq_unit_name)
        data = get_report(seq_unit_name)
        logger.debug(json.dumps(data, indent=3))
        mydict[seq_unit_name]["seq_proj_name"] = data["meta"]["seq_proj_name"]
        mydict[seq_unit_name]["seq_model"] = data["library_info"][0]["seq_model"]

        infos = ["library_name", "seq_unit_name", "raw_data_reads", "raw_data_size", "filter_data_reads", "filtered_data_size"]
        for info in infos:
            if info in data["filtered_fastq"]:
                mydict[seq_unit_name][info] = data["filtered_fastq"][info]

        nob_data = nob.Nob(data)
        n_stats = nob_data.find("stats")
        #print(n_stats)
        #print(n_stats[0])
        if len(n_stats) > 0:
            if nob_data[n_stats[0]][:] is not None:
                n_stats_value = nob_data[n_stats[0]][:]
            else:
                n_stats_value = -200
        else:
            n_stats_value = -100
        #print(n_stats_value)
        
        infos = [ "gc_divergence_coeff_r1_cg", "gc_divergence_coeff_r2_cg", "read_q20_read1", "read_q20_read2" ]
        for info in infos:
            mydict[seq_unit_name]["filtered_"+info] = "none"
            if "stats" in data["filtered_fastq"]:
                if data["filtered_fastq"]["stats"] is not None:
                    if info in data["filtered_fastq"]["stats"]:
                        mydict[seq_unit_name]["filtered_"+info] = data["filtered_fastq"]["stats"][info]

        mydict[seq_unit_name]["percent_filtered_reads"] = round(100 - float(mydict[seq_unit_name]["filter_data_reads"]) * 100 / float(mydict[seq_unit_name]["raw_data_reads"]), 2)
        mydict[seq_unit_name]["outputbases"] = int(data["filtered_fastq"]["stats"]["outputbases"])
        mydict[seq_unit_name]["outputbases[Gb]"] = round(float(mydict[seq_unit_name]["outputbases"]) / 1e9, 2)
        mydict[seq_unit_name]["complete_logical_amount (CLA)"] = round(float(mydict[seq_unit_name]["outputbases"]) / 1e6 / float(mydict[seq_unit_name]["genome-size-estimated-mb"]), 2)

        if "nt" in data["readqc"]["stats"]["sketch"]:
            full_tax = sorted( data["readqc"]["stats"]["sketch"]["nt"]["rows"], key= lambda d: int(d["matches"]), reverse=True)
            ###print(full_tax[0])
            ###full_tax_list = {k.split(":")[0]: k.split(":")[1] for k in full_tax[0]["taxonomy"].split(";")}
            full_tax_list = {}
            full_tax_array = full_tax[0]["taxonomy"].split(";")
            for ele in full_tax_array:
                k = ele.split(":")
                if len(k) > 1:
                    full_tax_list[k[0]] = k[1]
            #print(full_tax_list)
            mydict[seq_unit_name]["taxonomy species"] = full_tax_list["s"] if "s" in full_tax_list else "not_available"
            mydict[seq_unit_name]["taxonomy genus"]   = full_tax_list["g"] if "g" in full_tax_list else "not_available"
            mydict[seq_unit_name]["taxonomy family"]  = full_tax_list["f"] if "f" in full_tax_list else "not_available"
        else:
            mydict[seq_unit_name]["taxonomy species"] = "not_available"
            mydict[seq_unit_name]["taxonomy genus"] = "not_available"
            mydict[seq_unit_name]["taxonomy family"] = "not_available"

        infos = ["read_gc_mean", "read_gc_median", "read_gc_mode", "read_gc_std", "read_length_1", "read_length_2", "gc_divergence_coeff_r1_at", "gc_divergence_coeff_r1_atcg", "gc_divergence_coeff_r1_cg", "gc_divergence_coeff_r2_at", "gc_divergence_coeff_r2_atcg", "gc_divergence_coeff_r2_cg", "read_q20_read1", "read_q20_read2", "illumina_read_percent_contamination_phix", "illumina_read_percent_contamination_rna_spikein", "illumina_read_percent_contamination_dna_spikein", "illumina_read_percent_contamination_fosmid", "illumina_read_percent_contamination_microbes", "illumina_read_percent_contamination_rrna", "illumina_read_percent_contamination_plastid", "illumina_read_percent_contamination_mitochondrion", "illumina_read_percent_contamination_adapters", "illumina_read_percent_contamination_artifact_50bp", "illumina_read_percent_contamination_contaminants"]
        for info in infos:
            if info in data["readqc"]["stats"]:
                mydict[seq_unit_name][info] = data["readqc"]["stats"][info]

        # determine status for this library
        # all of the following needs to be true for the library to pass
        conditions = {}
        conditions["percent_filtered_reads"]    = float(mydict[seq_unit_name]["percent_filtered_reads"])             < configuration["ILLUMINA"].getfloat("percent_filtered_reads")
        conditions["complete_logical_amount"]   = float(mydict[seq_unit_name]["complete_logical_amount (CLA)"])      > configuration["ILLUMINA"].getfloat("complete_logical_amount")
        conditions["gc_divergence_coeff_r1_cg"] = float(mydict[seq_unit_name]["filtered_gc_divergence_coeff_r1_cg"]) < configuration["ILLUMINA"].getfloat("gc_divergence_coeff_r1_cg")
        conditions["gc_divergence_coeff_r2_cg"] = float(mydict[seq_unit_name]["filtered_gc_divergence_coeff_r2_cg"]) < configuration["ILLUMINA"].getfloat("gc_divergence_coeff_r2_cg")
        # since incorporation deteriorates toward the last cycles
        logger.info(json.dumps(conditions, indent=2))
        mydict[seq_unit_name]["status"] = "FAILED" if list(conditions.values()).count(False) > 0 else "PASSED"

        if mydict[seq_unit_name]["taxonomy genus"].lower() in mydict[seq_unit_name]["seq_proj_name"].lower():
            mydict[seq_unit_name]["warning"] = "none"
        else:
            mydict[seq_unit_name]["warning"] = "check taxonomy"

    # create pandas data fame
    dataframe = pandas.DataFrame.from_dict(mydict, orient="index")
    return dataframe


if __name__ == "__main__":
    # make sure environment variable LC_ALL is set, e.g. LC_ALL=en_US.UTF-8
    locale.setlocale(locale.LC_ALL, "en_US.UTF-8")

    # get logger
    logging.basicConfig(format="%(levelname)s:%(name)s:%(message)s", level="WARNING")
    logger = logging.getLogger(__name__)

    # get arguments
    args = get_arguments()

    # set log level and log output
    logger.setLevel(args.log.upper())
    # fh = logging.FileHandler('plantqc.log')
    # fh.setLevel(logging.DEBUG)
    # fh.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))

    # check arguments
    if args.output:
        if not args.output.endswith(".tsv"):
            logger.error(f"Output file name {args.output} is required to have suffix .tsv.")
            sys.exit(1)

    # get selection from configuration file
    config = configparser.ConfigParser()
    if not args.config or (args.config and not os.path.isfile(args.config)):
        # provide default values
        logger.warning("No configuration file provided, using default values!")
        config["ILLUMINA"] = {}
        config["ILLUMINA"]["percent_filtered_reads"] = '50'
        config["ILLUMINA"]["complete_logical_amount"] = '20.0'
        config["ILLUMINA"]["gc_divergence_coeff_r1_cg"] = '2.0'
        config["ILLUMINA"]["gc_divergence_coeff_r2_cg"] = '3.0'
    else:
        config_file = config.read(args.config)

    # print configuration
    if args.verbose:
        for section in config.sections():
            for item in config[section]:
                print(section, item, config[section][item])

    # get list of libraries to analyse
    if args.library_name:
        libraries = [args.library_name]
    elif args.input_file:
        if os.path.isfile(args.input_file):
            with open(args.input_file, "r", encoding="ascii") as file:
                libraries = [line.strip() for line in file]
                libraries = [line for line in libraries if line]  # remove blank lines
        else:
            sys.stderr.write("Error: file {args.input_file} does not exist. Exiting.\n")
            sys.exit(1)

    # calculate metrics
    df_list = []  # list of data frames
    platform_list = []
    for library in libraries:
        # get sequencing system
        linfo = get_library_info(library)
        platform = linfo["library_info"][0]["seq_platform"]
        platform_list.append(platform.lower())
        # get metrics
        if platform.lower() == "pacbio":
            dframe = pacbio_metrics(library, args.ccs)
        elif platform.lower() == "illumina":
            dframe = illumina_metrics(library, config)
        else:
            logger.error("unsupported platform %s", platform)
            sys.exit(1)
        df_list.append(dframe)

    # are all libraries from the same technology
    if platform_list.count(platform_list[0]) == len(platform_list):
        logger.info(f"same sequencing technology for all libraries: {platform_list[0]}")
    else:
        logger.error("not same sequencing technology for all libraries.")
        sys.exit(2)

    # write output
    df_total = pandas.concat(df_list)
    if args.verbose:
        print(df_total)

    if platform_list[0] == "pacbio":
        total_base_count = round(sum(df_total["sdm_raw_base_count"]) / 1e9, 1)  # in Gb
        total_read_count = sum(df_total["Subread Counts"])
        total_complete_logical_amount = sum(df_total["complete_logical_amount (CLA)"])
        summary = {}
        summary["Total number of reads"]  = locale.format_string('%d', total_read_count, grouping=args.grouping)
        summary["Total complete logical amount [X]"] = locale.format_string('%.1f', total_complete_logical_amount, grouping=args.grouping)
        summary["Total base count in Gb"] = str(total_base_count)
        summary["Target logical amount"] = max(df_total['target-logical-amount'])
        summary["Logical amount unit"] = df_total['logical-amount-units'][0]
        per_lib = df_total.groupby(["Library Name"])["sdm_raw_base_count"].agg(sum)
        per_lib = round(per_lib.astype(float) / 1e9, 1)
        df_per_lib = per_lib.to_frame()
        df_per_lib.columns = ["base count in Gb"]
        ds = pandas.DataFrame.from_dict(summary, orient='index')
        #ds.append(df_per_lib)
        print(ds)
   
    if args.grouping:
        numbers = list(df_total.select_dtypes(include=["int64"]).columns)
        for i in numbers:
            df_total[i] = df_total[i].apply(lambda x: locale.format_string("%d", x, grouping=True))
        reals = list(df_total.select_dtypes(include=["float64"]).columns)
        for i in reals:
            df_total[i] = df_total[i].apply(lambda x: locale.format_string("%.1f", x, grouping=True))

    if args.transpose:
        df_total = df_total.transpose()
        DO_INDEX = True
    else:
        DO_INDEX = False

    df_total.to_csv(args.output, sep="\t", index=DO_INDEX)

    if args.seq_unit:
        for i in df_total.columns:
            df = pandas.DataFrame(df_total[i])
            df.to_csv(args.output.replace("tsv", str(i) + ".tsv"), sep="\t", index=True)
            df.to_markdown(args.output.replace("tsv", str(i) + ".md"), index=True)

    if args.excel:
        df_total.to_excel(args.output.replace("tsv", "xlsx"), index=DO_INDEX)

    if args.markdown:
        # df_total.to_markdown(args.output.replace("tsv", "md"), index=DO_INDEX, tablefmt="grid")
        df_total.to_markdown(args.output.replace("tsv", "md"), index=DO_INDEX)

    if args.html:
        df_total.to_html(args.output.replace("tsv", "html"), index=DO_INDEX)

    if args.pdf:
        df_total.to_html(args.output.replace("tsv", "html"), index=DO_INDEX)
        pdfkit.from_file(args.output.replace("tsv", "html"), args.output.replace("tsv", "pdf"))

    logger.info(f"Created output {args.output}")

    if "jira_ticket" in args:
        try:
            jira = JIRA("https://issues.jgi.doe.gov", auth=(args.jira_user, args.jira_password))
        except JIRAError as e:
            logger.error(f"{e.status_code}: {e.text}")
        else:
            issue_comment = "{noformat}\n" + df_total.to_markdown() + "\n{noformat}\n"
            logger.info(issue_comment)
            comment_index = jira.add_comment(args.jira_ticket, issue_comment)

            issue = jira.issue(args.jira_ticket)
            trans_dict = {t["name"]: t["id"] for t in jira.transitions(issue)}
            logger.info(json.dumps(trans_dict))

            if "Complete" in trans_dict.keys():
                jira.transition_issue(issue, "Complete")
                logger.info("transition Complete.")
            else:
                logger.error("transition Complete not available.")
