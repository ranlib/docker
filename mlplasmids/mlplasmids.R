#!/usr/bin/Rscript

suppressMessages(library(mlplasmids))

usage =  "USAGE: Rscript run_mlplasmids.R ./path/to/assembly.fasta ./path/to/output.tab [prob_threshold] [species]"

args <- commandArgs(TRUE)

input_path <- args[1]
output_path <- args[2]
if(any(is.na(c(input_path, output_path)))){
    print(usage)
    stop()
}

# set the defaults
thresh <- ifelse(!is.na(args[3]), as.numeric(args[3]), .8)
species <- ifelse(!is.na(args[4]), args[4], "Escherichia coli")
print(paste("Threshold:", thresh))
print(paste("Species:", species))

example_prediction <- plasmid_classification(path_input_file = input_path,  prob_threshold=thresh, species = species)
if (is.null(example_prediction)){
    stop("Issue with mlplasmids; please try running interactively")
}

write.table(x=example_prediction, file=output_path, row.names=F, sep="\t")
