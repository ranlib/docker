#!/usr/bin/env Rscript

#args <- commandArgs(trailingOnly = TRUE)
#fn <- args[1]

## (1) set the input fasta file name. 
library(VirFinder)
inFaFile <- system.file("data", "contigs.fa", package="VirFinder")
#inFaFile <- system.file("data", fn, package="VirFinder")

## (2) prediction
predResult <- VF.pred(inFaFile)

#### (2.1) sort sequences by p-value in ascending order
predResult[order(predResult$pvalue),]

#### (2.2) estimate q-values (false discovery rates) based on p-values
predResult$qvalue <- VF.qvalue(predResult$pvalue)

#### (2.3) sort sequences by q-value in ascending order
predResult[order(predResult$qvalue),]

