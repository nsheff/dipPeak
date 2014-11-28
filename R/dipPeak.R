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
#' @references ' \url{http://github.com/sheffield}
## @import if you import any packages; here.
#' @docType package
#' @name dipPeak
#' @author Nathan Sheffield
NULL
################################################################################
#Now a few helper functions:
get_idxstats <- function(in_file) {
	if (Sys.which("samtools") == "") warning("samtools not installed; in get_idxstats());");
	stats = system(paste("samtools idxstats", in_file), intern=TRUE)
	stats_df = read.table(text=stats, sep="\t", col.names=c("contig", "length", "mapped", "unmapped"))
	#if (RESULT_STATS != 0 ) { cat("Warning: samtools failed; is samtools installed?\n"); }
}

#Takes a base path and creates a new, random directory for
#scratch writing. Use PREFIX to differentiate among software, if you want.
makescratchDir = function(scratchDirBase, PREFIX) {
	scratchDirBase = path.expand(scratchDirBase);
	#Make sure it has a trailing slash.
	if (substr(scratchDirBase,nchar(scratchDirBase), nchar(scratchDirBase)+1) != "/") {
		scratchDirBase = paste(scratchDirBase, "/", sep="");
	}
	repeat{#Find a unique folder
		randomFolder = floor(runif(1, min=0, max=1)*1e5); #so independent runs have their own private scratch space
		scratchDir = paste(scratchDirBase, PREFIX, randomFolder, "/",sep="") #include trailing slash!
		if(!file_test("-d", scratchDir))	{    	break;  	}
	}
	dir.create(scratchDir, showWarnings=FALSE, recursive=TRUE);
	return(scratchDir);
}

secToTime = function(timeInSec) {
	hr = timeInSec %/% 3600 #hours
	min = timeInSec %% 3600 %/% 60 #minutes
	sec = timeInSec %% 60 #seconds
	return(paste0(sprintf("%02d", hr), "h ", sprintf("%02d", min), "m ", sprintf("%02.01f", signif(sec, 3)), "s"))
}

tic <- function(gcFirst = TRUE, type=c("elapsed", "user.self", "sys.self")) {
	type <- match.arg(type)
	assign(".type", type, envir=baseenv())
	if(gcFirst) gc(FALSE)
	tic <- proc.time()[type]         
	assign(".tic", tic, envir=baseenv())
	invisible(tic)
	}

toc <- function() {
	type <- get(".type", envir=baseenv())
	toc <- proc.time()[type]
	tic <- get(".tic", envir=baseenv())
	timeInSec = as.numeric(toc-tic);
	message("<", secToTime(timeInSec), ">", appendLF=FALSE)
	invisible(toc)
}



################################################################################
#FUNCTION DOCUMENTATION - dipPeaks() main function
#' Estimate density and find peaks
#'
#' Calculates a smoothed track from input count data,
#' and identifies peaks in the data. Outputs the peaks and smoothed
#' track to file.
#'
#' @param bamFile	Input file in bam format
#' @param bigWigOut	Produce a smoothed output density file in bigwig format (using wigToBigWig)?
#' @param chromInfo		UCSC chromInfo.txt file (required only for bigWigOut)
#' @param scratchDirBase	Directory for scratch output; use local disk for optimal processing. dipPeak writes a significant amount of data to disk as it runs, to reduce memory use. This I/O takes time, and will be faster if you use a local disk or fast mount spot.
#' @param outDir	Final folder (the files to keep will be moved here from the scratch dir).
#' @param kernel		'e' for Epanechinov; 'g' for Gaussian; 'm' for Mexican Hat
#' @param cores		Add additional cores (up to 23) for faster processing. (Default:1)
#' @param perChromCutoff 	Boolean; If yes, calculates a per-chromosome; otherwise, defaults to genome-wide cutoff.
#' @param windowSize 	Choose the size of the smoothing window (default:50)
#' @param windowStep 	Choose the window step (5 will take a window every 5 bases) (default:5)
#' @param indexFile	Index file made from "samtools index" on the bam input file; deaults to bamFile.bai, or you can provide one.
#' @param retainTemp 	Do you want to retain temporary files?
#' @param limitChrom 	You can limit the smoothing to a list of chromosomes (defaults to everything in the bam file)
#' @export
#' @examples
#' dipPeaks("sequences.bam", bigWigOut="density.bw", outDir="outDir", cores=8);
dipPeaks = function(bamFile, bigWigOut=NULL, peakOut=NULL, chromInfo=NULL, scratchDirBase="~", outDir=NULL, kernel="e", cores=1, perChromCutoff=FALSE, windowSize=50, windowStep=5, indexFile=NULL, retainTemp=FALSE, limitChrom=NULL) {

### Process arguments
genomeWideCutoff = !perChromCutoff; 
useRle=TRUE; #RLE is now the only option.
#Check all options
if (!file.exists(bamFile)) { stop("Must define bamFile; file must exist."); }

if (!is.null(bigWigOut)) {
	if (is.null(chromInfo) || !file.exists(chromInfo)) { stop("To produce a bigwig file, you must define chromInfo (ucsc file for wigToBigWig); file must exist."); }
	if(Sys.which("wigToBigWig") == "") { stop("Can't find wigToBigWig; did you load the ucsc tools?") }
	#pass the tests, set this to true.
	outputDensity = TRUE;
} else {
	outputDensity = FALSE;
}

if(is.null(peakOut)) {
	if(is.null(bigWigOut)) {
		peakOut = "peaks"; #default file name
	} else {
		peakOut=paste0(bigWigOut, ".peaks");
	}
}

scratchDir=makescratchDir(scratchDirBase, PREFIX="dp");
if(file.access(scratchDir, mode=2) != 0) stop("Scratch dir [", scratchDir, "] not writable!\n");
on.exit(cleanUp(retainTemp, scratchDir))
if (is.null(outDir)) {
	outDir = getwd();
}
outDir = path.expand(outDir);
dir.create(outDir, showWarnings=FALSE, recursive=TRUE);
if(file.access(outDir, mode=2) != 0) stop("Output dir [", outDir, "] not writable!\n");


if (is.null(indexFile)) {
	indexFile=paste0(bamFile, ".bai"); #default index naming
}

### Kernel
#windowSize=50
x = (-(windowSize):(windowSize))/windowSize
ws=windowSize #bases on each side of the center base.
if(kernel == 'g') { #Gaussian
	kernelPrecise = dnorm(-(ws):(ws), mean=0, sd=ws/4) #kernelGaussian: a lower standard deviation will give you greater precision, higher will give a smoother density.
} else if (kernel == 'e') { #Epanechinokov
	kernelPrecise = (3/4) * (1-x^2) #kernelEpanechinokov
} else if (kernel == 'm') { #Mexican Hat
	
	sigma=1
	x = seq(from=-3.5, to=3.5, length.out=(ws*2+1))
	mh =  (1-((x^2)/(sigma^2)))*exp(-(x^2)/(2*sigma^2))
	mexicanHat = mh*ws/sum(abs(mh))
	mpe = .5 *(3/4) * (1-x^2) + .5 * mexicanHat
	kernelPrecise = mpe
} else {
	stop("Unrecognized kernel option.");
}
#You can also play with combining kernels.


#kernelNarrow = signif(mexicanHat, 4);
kernelWeights = signif(kernelPrecise,4)


### Output options for log files:
message("bamFile:", bamFile);
message("bigWigOut:", bigWigOut);
message("chromInfo:", chromInfo);
message("scratchDirBase:", scratchDirBase);
message("scratchDir:", scratchDir);
message("outDir:", outDir);
message("cores:", cores);
message("perChromCutoff:", perChromCutoff);
message("windowSize:", windowSize);
message("windowStep:", windowStep);
message("indexFile:", indexFile);
message("retainTemp:", retainTemp);
message("limitChrom:", limitChrom);
message("kernel:", kernel);



#tryCatch( { #for cleaning up temporary files.
cat("Copy input files to scratch space:", scratchDir, "...\n")
system(paste("cp", indexFile, scratchDir), wait=FALSE)
system(paste("cp", bamFile, scratchDir), wait=TRUE)
scratchBam = paste(scratchDir, basename(bamFile), sep="");

#genomeInfo = GRangesForUCSCGenome(genome="hg19")
idx = get_idxstats(scratchBam);
idx = idx[idx$mapped >0,]
genomeInfo = GRanges(seqnames=idx$contig, ranges=IRanges(start=rep(1, nrow(idx)), width=idx$length))

#chromlist = paste("chr", c(1:22, "X"), sep="")
#chromlist = paste("chr", c(20,21,22), sep="")
#Define which chromosomes to smooth; either from the bamFile or user input.
if(is.null(limitChrom)) {
	chromlist = as.list(as.vector(runValue(seqnames(genomeInfo))));
} else {
	if(is.list(limitChrom)) { chromlist=limitChrom; }
	else {chromlist = as.list(limitChrom); }
}

wiggleFiles = list();
tic();
#densityEstimate = list();
windowCount = list();
sampleFile = paste(scratchDir, "density.step.", windowSize, "-", windowStep, ".wig", sep="")
if(file.exists(sampleFile)) {
	system(paste("rm", sampleFile)) #clear room for sample file if it exists.
}
################################################################################
# Calculate density.
if (cores > 1) {
	library(parallel);
	windowCount = mclapply(chromlist, calculateDensity, genomeInfo,  scratchBam, windowSize, windowStep, kernelWeights, scratchDir, useRle, outputDensity, genomeWideCutoff, mc.cores=cores)
} else {
	windowCount = lapply(chromlist, calculateDensity,  genomeInfo,  scratchBam, windowSize, windowStep, kernelWeights, scratchDir, useRle, outputDensity, genomeWideCutoff)
}

names(windowCount) = chromlist

if (genomeWideCutoff) {
	cat("Calculating genome-wide cutoff...");
	#t = scan("/scratch/NCS/chr1.density.b.50-5.wig", skip=1, allowEscapes=TRUE)
	#Slurp in the sampled file:
	system(paste("cat `ls ", scratchDir,"*.sample*` > ", scratchDir, "/densitySample.wig", sep=""))
	sampleScores = scan(file=paste(scratchDir, "/densitySample.wig", sep=""), allowEscapes=TRUE)
	#the step ones are small enough that this can probably be fast.
	looseQuantile = fitGammaCutoff(sampleScores, 0.9);
	if (is.nan(looseQuantile)) { stop("looseQuantile is NaN; try a different kernel"); }
	cat("Cutoff:",looseQuantile,"\n");
	rm(sampleScores); gc(); toc();		#free up some memory
	cat("Finding peaks...");

	if (cores > 1) {
		result = mclapply(chromlist, defineAllPeaks, scratchDir, windowSize, windowStep, windowCount, looseQuantile, mc.cores=cores)
	} else {
		result = lapply(chromlist, defineAllPeaks, scratchDir, windowSize, windowStep, windowCount, looseQuantile)
	}
}
peakOutDips = paste0(peakOut, ".dips.bed");
peakOutBroad = paste0(peakOut, ".broad.bed");
message("\nPeak finding complete.")
system(paste("cat `ls ", scratchDir,"chr*.peaks.dips.bed` > ", outDir, "/", peakOutDips, sep=""))
system(paste("cat `ls ", scratchDir,"chr*.peaks.b.bed` > ", outDir, "/", peakOutBroad, sep=""))
BEDSORT_RESULT = system(paste("bedSort ", outDir, "/", peakOutDips, " ", outDir, "/", peakOutDips, sep=""));
BEDSORT_RESULT2 = system(paste("bedSort ", outDir, "/", peakOutBroad, " ", outDir, "/", peakOutBroad, sep=""));

if (BEDSORT_RESULT + BEDSORT_RESULT2 > 0) {
	warning("Warning; Bedsort failed. Is bedsort installed?\n");
}


outputDensity_RESULT=0;
if (!is.null(bigWigOut)) {
	#makeBigWig = paste("cat ", paste(wiggleFiles, collapse=" "), " | wigToBigWig -clip stdin chromInfo.hg19.txt ", scratchOutfile, "; mv ", scratchOutfile, outDir)scratchDir
scratchOutfile = paste(scratchDir, basename(bigWigOut), sep="");
	makeBigWig = paste("cat `ls ", scratchDir,"chr*.dens*.wig` | wigToBigWig -clip stdin ", chromInfo, " ", scratchOutfile, "; mv ", scratchOutfile, " ", outDir, sep="")
	cat("Building density bigwig... Command:", makeBigWig, "\n");
	tic();outputDensity_RESULT = system(makeBigWig, wait=TRUE);toc();
}
if (!retainTemp & outputDensity_RESULT == 0) { #success bigwig!
	message("\nDipPeak run finished, output folder [", outDir, "]")
} else {
	message("\nbigWigOut process failed! Couldn't create a bigwig file with wigToBigWig.");
}
### closing tryCatch code.
#}, interrupt= function(x) { print(x) },
#error = function(x) { print(x) },
#finally = cleanUp(retainTemp, scratchDir)) #end tryCatch
}#END PARENT FUNCTION

################################################################################
#FUNCTION DOCUMENTATION - dipPeaks() main function
#' Just a display of what the different kernels look like.
#' @param windowSize	window size on the two sides of the center.
#@export
plotKernels = function(windowSize=50) {
	ws=windowSize #bases on each side of the center base.
	x = (-(ws):(ws))/ws
	kernelGaussian = dnorm(-(ws):(ws), mean=0, sd=ws/4) #kernelGaussian
	kernelEpanechinokov = (3/4) * (1-x^2) #kernelEpanechinokov
	sigma=1
	x = seq(from=-3.5, to=3.5, length.out=(ws*2+1))
	mh =  (1-((x^2)/(sigma^2)))*exp(-(x^2)/(2*sigma^2))
	kernelMexicanHat = mh*ws/sum(abs(mh))
	mpe = .5 *kernelEpanechinokov + .5 * kernelMexicanHat
	par(pty="s");
	ymin = min(kernelGaussian, kernelEpanechinokov, kernelMexicanHat, mpe)
	plot(kernelMexicanHat/max(kernelMexicanHat), type="l", main="Kernels", ylab="Signal Density", xlab="Nucleotides", lwd=2, ylim=c(ymin,1))
	lines(kernelGaussian/max(kernelGaussian), col="blue", lwd=2)
	lines(kernelEpanechinokov/max(kernelEpanechinokov), col="red", lwd=2)
	lines(mpe/max(mpe), col="maroon", lwd=2)
	legend('topright', c("mexican hat", "gaussian", "epanechinokov", "epa + mh"), col=c("black", "blue", "red", "maroon"), lwd=2, cex=.7)
}
#plotKernels()

#' Internal function to remove temporary files on interrupt or completion.
cleanUp = function(retainTemp, scratchDir) {
	if (!retainTemp) {
		message('\nCleaning up temporary files in scratchDir:\t', scratchDir, '...');
		unlink(scratchDir, recursive=TRUE);
	} else {
		message('\nTemporary files have been retained in scratchDir:\t', scratchDir);
	}
}












