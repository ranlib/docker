#!/usr/bin/env python
"""
get pacbio plant QC report
adapted from
Vasanth Singan
"""
import sys
import os
import argparse
from operator import itemgetter
from collections import OrderedDict
import csv
import requests
import jgidb

desc = """PBstats.py - takes list of libraries and creates table of PacBio stats"""

def get_arguments():
    """get arguments"""
    parser = argparse.ArgumentParser(description=desc, prog=os.path.basename(__file__), formatter_class=lambda prog: argparse.HelpFormatter(prog, max_help_position=80))
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--file", "-f", dest="input_file", type=str, help="text file with library names, one per line")
    group.add_argument("--lib", "-l", dest="library_name", type=str, help="name of library")
    parser.add_argument("--output", "-o", dest="output", type=str, help="output filename (default: %(default)s)", default="Lib_QC.txt")
    parser.add_argument("--verbose", "-v", dest="verbose", help="modify output verbosity (default: %(default)s)", action="store_true")
    return parser.parse_args()


def top_blast_hit(parsed_blastout: str) -> str:
    """get top blast hit"""
    family = "None"
    if os.path.isfile(parsed_blastout):
        with open(parsed_blastout, "r") as f:
            lines = [line.split() for line in f]
        lines.sort(key=itemgetter(2), reverse=True)
        firstline = lines[0]
        taxonomy = firstline[-1].strip().split(";")
        if len(taxonomy) > 4:
            family = taxonomy[4]
        else:
            family = taxonomy[0]
    return family


def top_sketch_hit(sketch_file: str) -> str:
    """get top sketch hit"""
    ln = 0
    family = "None"
    if os.path.isfile(sketch_file):
        with open(sketch_file, "r") as f:
            for line in f:
                ln += 1
                if ln == 4 and len(line) > 2:
                    try:
                        # print line
                        family = line.split(";")[-3].replace("f:", "")
                    except IndexError:
                        family = "None"
    return family


def pb_metrics(libs: [], output: str, verbose: bool):
    """query library in database"""
    sqlLibs = """
SELECT 
   l.library_name as 'Library Name', 
   s.seq_unit_name as 'Seq Unit Name', 
   l.tax_family as 'Sample Family', 
   s.raw_reads_count as 'Subread Counts', 
   ps.p0 as P0, ps.p1 as P1, ps.p2 as P2,
   ps.pf_read_count as 'PostFilter Read Count',
   ps.insert_read_length as 'Reads of Insert Length',
  (select concat(rpf.fs_location,rpf.file_name) from rqc_pipeline_files rpf where rpf.rqc_pipeline_queue_id = rpq.rqc_pipeline_queue_id
  and rpf.is_production = 1
  and rpf.file_type = 'pacbioSubreadQc.stats.readQual.Longest_Subread.hist' limit 0,1) as File1,
  (select concat(rpf.fs_location,rpf.file_name) from rqc_pipeline_files rpf where rpf.rqc_pipeline_queue_id = rpq.rqc_pipeline_queue_id and rpf.is_production = 1
  and rpf.file_type = 'sketch nt' limit 0,1) as File2
  from seq_units s
  inner join library_info l on s.rqc_library_id = l.library_id
  left join pacbio_stats ps on s.seq_unit_id = ps.seq_unit_id
  inner join rqc_pipeline_queue rpq on s.seq_unit_id = rpq.seq_unit_id and rpq.rqc_pipeline_type_id = 49 and rpq.rqc_pipeline_status_id in (14,31)
  where l.library_name IN ("""

    i = 0
    paramsLibs = []
    for lib in libs:
        i += 1
        paramsLibs.append(lib)
        sqlLibs = sqlLibs + "'" + lib + "', "

    sqlLibs = sqlLibs[:-2]
    sqlLibs = sqlLibs + ")"

    mydict = OrderedDict()

    # get information from database via REST
    mykeys = ["Subread Counts", "P2", "Library Name", "P1", "Reads of Insert Length", "P0", "PostFilter Read Count", "File1", "Seq Unit Name", "Sample Family", "File2" ]
    for lib in libs:
        report = "https://rqc.jgi-psf.org/api/rqcws/pbstats/" + lib
        #print(report)
        try:
            response = requests.get(report)
        except requests.exceptions.RequestException as e:
            print(e)
            continue
        data = response.json()
        for seq_unit in data["seq_units"]:
            name = seq_unit["Seq Unit Name"]
            #for key in mykeys:
            #    mydict[name][key] = seq_unit[key]
        
    # get information from database
    for row in jgidb.queryDb(dbName="RQC", sql=sqlLibs):  # , sqlParams=paramsLibs):
        if "subread" in row["Seq Unit Name"] or "ccs" in row["Seq Unit Name"]:
            if row["Seq Unit Name"] not in mydict:
                mydict[row["Seq Unit Name"]] = row
        else:
            print("<W> pb_metrics: neither subreads nor ccs reads available for library ")

    # calculate some derived quantities
    for key,value in mydict.items(): # key = seq_unit
        Norm = value["P0"] + value["P1"] + value["P2"]
        if Norm > 0:
            value["P0[%]"] = round( value["P0"] * 100 / Norm, 2)
            value["P1[%]"] = round( value["P1"] * 100 / Norm, 2)
            value["P2[%]"] = round( value["P2"] * 100 / Norm, 2)
        else:
            value["P0[%]"] = value["P0"]
            value["P1[%]"] = value["P1"]
            value["P2[%]"] = value["P2"]

    # check out family
    for key,value in mydict.items(): # key = seq_unit
        mydict[key]["Family"] = top_sketch_hit(value["File2"])
        # family = top_blast_hit(value['File2'])
        if value["Sample Family"] == "None":
            mydict[key]["Correct Family"] = "No Info"
        else:
            mydict[key]["Correct Family"] = "Yes" if value["Sample Family"] == value['Family'] else "No"

    # get raw base count
    for lib in libs:
        report = 'https://rqc.jgi-psf.org/api/fullreport/report/' + lib
        #print(report)
        try:
            response = requests.get(report)
        except requests.exceptions.RequestException as e:
            print(e)
            continue
        data = response.json()
        for library in data["library_info"]:
            seq_unit = library["seq_unit_name"]
            mydict[seq_unit]["sdm_raw_base_count"] = library["sdm_raw_base_count"]
            mydict[seq_unit]["Raw base count [Gb]"] = round( float(library["sdm_raw_base_count"]) / 1e9, 1 )


    # debugging print out
    if verbose:
        print("<I> pb_metrics: mydict = ")
        print(mydict)

    # write output
    invalid = [ "P0", "P1", "P2", "File1", "File2" ]
    header = list(next(iter(mydict.items()))[1].keys())
    #header = list(list(mydict.values())[0].keys())
    new_header = [ x for x in header if x not in invalid ]
    with open(output, "w") as out:
        writer = csv.DictWriter(out, fieldnames=new_header, delimiter='|')
        writer.writeheader()
        for key,value in mydict.items(): # key = seq_unit
            record = { i:j for i,j  in value.items() if i not in invalid }
            record["Subread Counts"] = '{:,}'.format(record["Subread Counts"])
            record["PostFilter Read Count"] = '{:,}'.format(record["PostFilter Read Count"])
            record["Reads of Insert Length"] = '{:,}'.format(int(record["Reads of Insert Length"]))
            record["sdm_raw_base_count"] = '{:,}'.format(int(record["sdm_raw_base_count"]))
            record["Raw base count [Gb]"] = '{:,}'.format(record["Raw base count [Gb]"])
            writer.writerow(record)

    print("Created output " + output)


if __name__ == "__main__":
    args = get_arguments()
    if args.library_name:
        libs = [ args.library_name ]
        pb_metrics(libs, args.output, args.verbose)
    elif args.input_file:
        if os.path.isfile(args.input_file):
            with open(args.input_file, "r") as file:
                libs = [ line.strip() for line in file ]
                pb_metrics(libs, args.output, args.verbose)
        else:
            sys.stderr.write("Error: libs file '%s' does not exist. Exiting.\n" % args.input_file)
            sys.exit(1)
