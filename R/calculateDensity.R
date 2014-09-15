
calculateDensity = function(chrom, genomeInfo, scratchBam, windowSize, windowStep, kernelWeights,  SCRATCHDIR, USE_RLE, OUTPUT_DENSITY, GENOME_WIDE_CUTOFF) {
	cat("\nChrom:", chrom, "\t");
	cat("Loading sequences...\t");
	wholeChrom =  genomeInfo[seqnames(genomeInfo) == chrom]
	what = c("rname", "strand", "pos")
	param = ScanBamParam(which = wholeChrom, what=what)
	sb = scanBam(scratchBam, param=param);toc();
	cat("Calculating density...\t");
	name <- names(sb)
	#readsGR<-GRanges(seqnames=as.vector(sb[[name]]$rname),IRanges(start=sb[[name]]$pos,width=1),strand=sb[[name]]$strand);
	#wholeChrom = GRangesForUCSCGenome(genomeBuild, chrom=as.character(runValue(seqnames(singleChromGR))))
	#wholeChromData <- import.bw(bamFile, BigWigSelection(wholeChrom), asRangedData=F);
	#cov = coverage(readsGR) #resize(readsGR,1)

	#using IRanges:
	reads = IRanges(start=sb[[name]]$pos,width=1)
	cov = coverage(reads);

	if (USE_RLE) { #if using GenomicRanges, sub "cov" with "cov[[1]]" to access the innards. If using IRanges, use "cov".
		completeViews = successiveViews(cov, width=rep((windowSize*2+1),length(cov)/windowStep), gapwidth=-((windowSize*2+1)-windowStep))
		windowMulsBroad = viewMuls(completeViews, na.rm=FALSE, kernelWeights);toc();
		#object.size(completeViewsR)/1000/1000
	} else {
		completeViews = successiveViews(as.vector(cov), width=rep((windowSize*2+1),length(cov)/windowStep), gapwidth=-((windowSize*2+1)-windowStep))
		windowMulsBroad = viewMuls(completeViews, na.rm=FALSE, kernelWeights);toc();
		#object.size(completeViews)/1000/1000
	}
	windowCount = length(completeViews);
	rm(completeViews); gc();
	densityEstimate = signif(windowMulsBroad, digits=3); toc();
	rm(windowMulsBroad); gc();
	if (OUTPUT_DENSITY || GENOME_WIDE_CUTOFF) {
		#If the user requests to save the the density to disk, save here. Otherwise they will be discarded.
		#If the user requests a genome-wide cutoff, even without wanting to save the output, we still have
		#to put them onto the disk temporarily, until we can calculate the genome-wide cutoff; then we will
		#read them in again, and can delete the temp files.
		cat("Printing density...\t");
		wiggleFileName = paste(SCRATCHDIR, chrom, ".density.", windowSize, "-", windowStep, ".wig", sep="")
		write(paste("fixedStep  chrom=", chrom, " start=",windowSize+1,"  step=", windowStep, "  span=", 1, sep=""), file=wiggleFileName, sep="\t")
		#tic(); write.table(densityEstimates2, row.names=FALSE, col.names=FALSE, quote=FALSE, file=wiggleFileName, append=TRUE); toc();
		cwriteR(as.vector(densityEstimate), as.numeric(length(densityEstimate)), file=wiggleFileName); 

toc();
	}
	if (GENOME_WIDE_CUTOFF) {
		#If the user requests a genome-wide cutoff, then we also will spit out a sample (default: 10%) 
		#of every chromosome into a temporary file. This is faster than using all the data, 
		#and more accurate than using only a single chromosome, which may be biased.
		cat("Storing sample of density for genome-wide cutoff calculation...\t");
		sampleFileName = paste(SCRATCHDIR, chrom, ".sample.", windowSize, "-", windowStep, ".wig", sep="")
		file.create(sampleFileName);
		cwriteStepR(as.vector(densityEstimate), as.numeric(length(densityEstimate)), file=sampleFileName, step=10); toc();
	}

	if (!GENOME_WIDE_CUTOFF) {
		#The alternative is to use a different cutoff for each chromosome. This can lead to problems if
		#some chromosomes have strange distributions or less data.
		looseQuantile = quantile(densityEstimate[densityEstimate > 0], probs=c(.95))
		windowCoords = seq(from=51, by=5, length.out=windowCount)
		peakCoords = defineBroadPeaks(chrom, densityEstimate, windowCoords, SCRATCHDIR, looseQuantile,windowStep);
		defineNarrowPeaks(chrom, densityEstimate, windowCoords, SCRATCHDIR, peakCoords);
	}
	return(windowCount);
}
