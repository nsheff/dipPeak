#PACKAGE DOCUMENTATION
#' Parzen density estimation plus peak divider based on dips
#'
#' dipPeak identifies peaks in high-throughput sequencing data.
#' It is designed for use with DNase-seq data, but could be useful for
#' other data types as well. It does two tasks: First, identifies broad
#' peaks (hypersensitive regions), and then it divides these peaks into
#' smaller regions of interest to be considered separately. 
#'
#' It also has the advantage of allowing user control over the resolution
#' of the output density. This can allow the user to determine the tradeoff
#' between disk space and resolution.
#' 
#' @references This has not been published but check out my web page
#' \url{http://www.nathansheffield.com}
## @import if you import any packages; here.
#' @docType package
#' @name dipPeak
#' @author Nathan Sheffield
NULL

#I made it so you can do this stuff on RLE vectors, 
#which cuts down the memory use to about 1/3 in my tests,
#and speeds it up, too.
#USE_RLE = TRUE;

#At the cost of more CPU use and some scratch space, you can estimate the cutoff
#using all data, instead of a different cutoff based on each
#individual chromosome. This is highly recommended, since sometimes single 
#chromosome cutoffs can lead to weird results.
#GENOME_WIDE_CUTOFF = TRUE;

#wiggle density files can be big. If you don't want to keep them
#you don't have to.
#OUTPUT_DENSITY = TRUE;
tic = function() {}
toc = function() {}


dipPeakOptions = data.frame(
BAMFILE="encodeRegions.bam",
INDEXFILE="encodeRegions.bam.bai",
SCRATCHDIRBASE="scratch",
BIGWIGOUT="", 
OUTDIR="finalOutput", 
CORES=1,
windowSize=50,
windowStep=5,
USE_RLE=TRUE,
OUTPUT_DENSITY=TRUE,
GENOME_WIDE_CUTOFF=TRUE);



get_idxstats <- function(in_file) {
	stats = system(paste("samtools idxstats", in_file), intern=TRUE)
	stats_df = read.table(text=stats, sep="\t", col.names=c("contig", "length", "mapped", "unmapped"))
	#if (RESULT_STATS != 0 ) { cat("Warning: samtools failed; is samtools installed?\n"); }
}


#################################
#Takes a base path and creates a new, random directory for
#scratch writing. Use PREFIX to differentiate among software, if you want.
makeScratchDir = function(SCRATCHDIRBASE, PREFIX) {
	SCRATCHDIRBASE = path.expand(SCRATCHDIRBASE);
	#Make sure it has a trailing slash.
	if (substr(SCRATCHDIRBASE,nchar(SCRATCHDIRBASE), nchar(SCRATCHDIRBASE)+1) != "/") {
		SCRATCHDIRBASE = paste(SCRATCHDIRBASE, "/", sep="");
	}
	repeat{#Find a unique folder
		randomFolder = floor(runif(1, min=0, max=1)*1e5); #so independent runs have their own private scratch space
		SCRATCHDIR = paste(SCRATCHDIRBASE, PREFIX, randomFolder, "/",sep="") #include trailing slash!
		if(!file_test("-d", SCRATCHDIR))	{    	break;  	}
	}
	dir.create(SCRATCHDIR, showWarnings=FALSE, recursive=TRUE);
	return(SCRATCHDIR);
}
################################


#FUNCTION DOCUMENTATION
#' Estimate density and find peaks
#'
#' Calculates a smoothed track from input count data,
#' and identifies peaks in the data. Outputs the peaks and smoothed
#' track to file.
#'
#'
#'
#' @param BAMFILE	Input file in bam format
#' @param INDEXFILE	Index file made from "samtools index" on the bam input file
#' @param SCRATCHDIRBASE	Directory for scratch output; use local disk for optimal processing. dipPeak writes a significant amount of data to disk as it runs, to reduce memory use. This I/O takes time, and will be faster if you use a local disk or fast mount spot.
#' @param BIGWIGOUT	Smoothed output density file in bigwig format (if you request an output density).
#' @param OUTDIR	Final folder (the files to keep will be moved here from the scratch dir).
#' @param CORES		Add additional cores (up to 23) for faster processing. (Default:1)
#' @param USE_RLE	Boolean; offers speed and memory improvement
#' @param GENOME_WIDE_CUTOFF 	Boolean; If yes, calculates a genome-wide cutoff for peak detection. Otherwise, it does this on a per-chromosome basis
#' @param OUTPUT_DENSITY 	Boolean; If yes, produces a bigwig density track.
#' @param windowSize 	Choose the size of the smoothing window (default:50)
#' @param windowStep 	Choose the window step (5 will take a window every 5 bases) (default:5)
#' @param CHROMINFO		UCSC chromInfo.txt file
#' @param RETAIN_TEMP 	Do you want to retain temporary files?
#' @param LIMIT_CHROM 	You can limit the smoothing to a list of chromosomes (defaults to everything in the bam file)
#' @export
#' @examples
#' dipPeaks("sequences.bam", "sequences.bam.bai", "scratch", "density.bw", "OUTDIR", CORES=8);
#' dipPeaks("encodeRegions.bam", "encodeRegions.bam.bai", "scratch", "bw.bw", "finalOutput", 8)
dipPeaks = function(BAMFILE, INDEXFILE, SCRATCHDIRBASE, BIGWIGOUT, OUTDIR, CORES=1, USE_RLE=TRUE, GENOME_WIDE_CUTOFF=TRUE, OUTPUT_DENSITY=TRUE, windowSize=50, windowStep=5, CHROMINFO, RETAIN_TEMP=FALSE, LIMIT_CHROM=NULL) {

#Check all options
if (OUTPUT_DENSITY) {
	if (!exists("CHROMINFO")) || !file.exists(CHROMINFO)) { stop("Must define CHROMINFO; file must exist."); }
}

SCRATCHDIR=makeScratchDir(SCRATCHDIRBASE, PREFIX="dp");
if(file.access(SCRATCHDIR, mode=2) != 0) stop("SCRATCHDIR [", SCRATCHDIR, "] not writable!\n");

OUTDIR = path.expand(OUTDIR);
dir.create(OUTDIR, showWarnings=FALSE, recursive=TRUE);
if(file.access(OUTDIR, mode=2) != 0) stop("OUTDIR [", OUTDIR, "] not writable!\n");

cat("Copy input files to scratch space:", SCRATCHDIR, "...\n")
system(paste("cp", INDEXFILE, SCRATCHDIR), wait=FALSE)
system(paste("cp", BAMFILE, SCRATCHDIR), wait=TRUE)
scratchBam = paste(SCRATCHDIR, basename(BAMFILE), sep="");
scratchOutfile = paste(SCRATCHDIR, basename(BIGWIGOUT), sep="");


#h=150 #bases on each side of the center base.
#sigma=1
#x = seq(from=-3.5, to=3.5, length.out=(h*2+1))
#mh =  (1-((x^2)/(sigma^2)))*exp(-(x^2)/(2*sigma^2))
#mexicanHat = mh*h/sum(abs(mh))

#windowSize=50
x = (-(windowSize):(windowSize))/windowSize
precomputedEpanechinokov = (3/4) * (1-x^2)
#precomputedGaussian = dnorm(-(h):(h), mean=0, sd=h/4) #a lower standard deviation will give you greater precision, higher will give a smoother density.
#mpe = .5 *precomputedEpanechinokov + .5 * mexicanHat

#windowStep=5
#kernelNarrow = signif(mexicanHat, 4);
kernelWeights = signif(precomputedEpanechinokov,4)


#genomeInfo = GRangesForUCSCGenome(genome="hg19")
idx = get_idxstats(scratchBam);
idx = idx[idx$mapped >0,]
genomeInfo = GRanges(seqnames=idx$contig, ranges=IRanges(start=rep(1, nrow(idx)), width=idx$length))


#chromlist = paste("chr", c(1:22, "X"), sep="")
#chromlist = paste("chr", c(20,21,22), sep="")
#Define which chromosomes to smooth; either from the bamfile or user input.
if(is.null(LIMIT_CHROM)) {
	chromlist = as.list(as.vector(runValue(seqnames(genomeInfo))));
} else {
	if(is.list(LIMIT_CHROM)) { chromlist=LIMIT_CHROM; }
	else {chromlist = as.list(LIMIT_CHROM); }
}

wiggleFiles = list();
tic();
#densityEstimate = list();
windowCount = list();
sampleFile = paste(SCRATCHDIR, "density.step.", windowSize, "-", windowStep, ".wig", sep="")
system(paste("rm", sampleFile))

if (CORES > 1) {
	library(multicore);
	windowCount = mclapply(chromlist, calculateDensity, genomeInfo,  scratchBam, windowSize, windowStep, kernelWeights, SCRATCHDIR, USE_RLE, OUTPUT_DENSITY, GENOME_WIDE_CUTOFF, mc.cores=CORES)
} else {
	windowCount = lapply(chromlist, calculateDensity,  genomeInfo,  scratchBam, windowSize, windowStep, kernelWeights, SCRATCHDIR, USE_RLE, OUTPUT_DENSITY, GENOME_WIDE_CUTOFF)
}

names(windowCount) = chromlist

if (GENOME_WIDE_CUTOFF) {
	cat("Calculating genome-wide cutoff...");
	#t = scan("/scratch/NCS/chr1.density.b.50-5.wig", skip=1, allowEscapes=TRUE)
	#Slurp in the sampled file:
	system(paste("cat `ls ", SCRATCHDIR,"*.sample*` > ", SCRATCHDIR, "/densitySample.wig", sep=""))
	sampleScores = scan(file=paste(SCRATCHDIR, "/densitySample.wig", sep=""), allowEscapes=TRUE)
	#the step ones are small enough that this can probably be fast.
	looseQuantile = fitGammaCutoff(sampleScores, 0.9);
	cat("Cutoff:",looseQuantile,"\n");
	rm(sampleScores); gc(); toc();		#free up some memory
	cat("Finding peaks...");

	if (CORES > 1) {
		result = mclapply(chromlist, defineAllPeaks, SCRATCHDIR, windowSize, windowStep, windowCount, looseQuantile, mc.cores=CORES)
	} else {
		result = lapply(chromlist, defineAllPeaks, SCRATCHDIR, windowSize, windowStep, windowCount, looseQuantile)
	}
}

system(paste("cat `ls ", SCRATCHDIR,"chr*.peaks.dips.bed` > ", OUTDIR, "/peaks.dips.bed", sep=""))
system(paste("cat `ls ", SCRATCHDIR,"chr*.peaks.b.bed` > ", OUTDIR, "/peaks.b.bed", sep=""))
BEDSORT_RESULT = system(paste("bedSort ", OUTDIR, "/peaks.dips.bed ", OUTDIR, "/peaks.dips.bed", sep=""));
BEDSORT_RESULT2 = system(paste("bedSort ", OUTDIR, "/peaks.b.bed ", OUTDIR, "/peaks.b.bed", sep=""));

if (BEDSORT_RESULT + BEDSORT_RESULT2 > 0) {
	cat("Warning; Bedsort failed. Is bedsort installed?\n");
	}


toc();
if (OUTPUT_DENSITY) {
	#makeBigWig = paste("cat ", paste(wiggleFiles, collapse=" "), " | wigToBigWig -clip stdin chromInfo.hg19.txt ", scratchOutfile, "; mv ", scratchOutfile, OUTDIR)SCRATCHDIR
	makeBigWig = paste("cat `ls ", SCRATCHDIR,"chr*.dens*.wig` | wigToBigWig -clip stdin ", CHROMINFO, " ", scratchOutfile, "; mv ", scratchOutfile, " ", OUTDIR, sep="")
	cat(makeBigWig);
	cat("Building density bigwig...\n");
	tic();OUTPUT_DENSITY_RESULT = system(makeBigWig, wait=TRUE);toc();
}

if (!RETAIN_TEMP & OUTPUT_DENSITY_RESULT == 0) { #success bigwig!
	system(paste("rm -rf", SCRATCHDIR))
}
}#END PARENT FUNCTION

