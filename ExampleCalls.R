#!/usr/bin/Rscript
#############################################
# Parzen density estimator - DipSeq
# By Nathan Sheffield, Duke University, 2013
#############################################
# This script reads a BAM file for a high-throughput sequencing experiment, and then 
# uses a density estimator (in C++) to calculate a coverage density, outputs to 
# wiggle (then converted to bigwig), and also draws a cutoff to define peaks in signal.
#
# It relies on a modified IRanges function I wrote.
# It also uses my R package for fast C output called cwrite.
#
# install.packages("~/IRanges", repos=NULL);
#
usage=' dipSeq.R <reads.bam> <outFolder> <output.bw> <scratchDir>';
##############################################
#genomeBuild="hg19"
##############################################
BAMFILE = "sequence.final.bam";
OUTFOLDER = "out/";
BIGWIGOUT = "dipSeq.bw";
SCRATCHDIRBASE = "~/";
CORES = 4;
#Command-line reading.
##########################################
source("~/bin/funcCommon.R")
#The command-line args are listed here in order:
variables = c(
"BAMFILE", 
"OUTFOLDER",
"BIGWIGOUT",
"SCRATCHDIRBASE",
"CORES"
)
#variables = data.frame(BAMFILE="sequence.final.bam", BIGWIGOUT="dipSeq.bw", SCRATCHDIRBASE="/scratch/dipSeq", GENOME="hg19")

ARGS = loadArgs(variables, minimum_args=length(variables)-3);
##########################################
#RDATA.DIR="rData/";
#FIGURE.DIR = "fig/";
#dir.create(FIGURE.DIR, showWarnings=FALSE);
#dir.create(RDATA.DIR, showWarnings=FALSE);

##############################################
#Load up the libraries
library(GenomicRanges)
library(rtracklayer)
library(Rsamtools)
library(cwrite)
library(IRangesMultiplication)
#OR:
library(dipPeak);


source("~/bin/funcPeaks.R")
#library(cwrite, lib.loc="~/cwrite");


if (ARGS$CORES > 1) {
	library(multicore);
}

##############################################
indexFile = paste(ARGS$BAMFILE, ".bai", sep="")
if(!file_test("-f", indexFile))	{    	
	system(paste("samtools index ", ARGS$BAMFILE, sep=""));
}

outFolder = ARGS$OUTFOLDER;

dipPeak(BAMFILE, indexFile, OUTFOLDER, CORES, USE_RLE=TRUE, GENOME_WIDE_CUTOFF=TRUE, OUTPUT_DENSITY=TRUE) 

source("Rncs/dipPeak/R/funcPeaks.R")
source("Rncs/dipPeak/R/definePeaks.R")
source("Rncs/dipPeak/R/dipPeak.R")
source("Rncs/dipPeak/R/calculateDensity.R")

library(dipPeak)
BAMFILE = "encodeRegions.bam";
INDEXFILE = "encodeRegions.bam.bai";
BAMFILE = "sequence.final.bam";
INDEXFILE = "sequence.final.bam.bai";
SCRATCHDIRBASE = "scratch";
BIGWIGOUT = "dipPeak.bw";
OUTDIR = "finalOutput"; #no trailing slash
CORES = 8;
USE_RLE=TRUE
GENOME_WIDE_CUTOFF=TRUE
OUTPUT_DENSITY=TRUE
windowSize=50
windowStep=5

dipPeak(BAMFILE, INDEXFILE, SCRATCHDIRBASE, BIGWIGOUT, OUTFOLDER, CORES, USE_RLE=TRUE, GENOME_WIDE_CUTOFF=TRUE, OUTPUT_DENSITY=TRUE, windowSize=50, windowStep=5)




