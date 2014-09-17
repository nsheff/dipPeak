
calculateDensity = function(chrom, genomeInfo, scratchBam, windowSize, windowStep, kernelWeights,  scratchDir, useRle, outputDensity, genomeWideCutoff) {
	message("\n[", chrom, "]\t", appendLF=FALSE);
	message("Loading sequences...\t", appendLF=FALSE);
	wholeChrom =  genomeInfo[seqnames(genomeInfo) == chrom]
	what = c("rname", "strand", "pos")
	param = ScanBamParam(which = wholeChrom, what=what)
	sb = scanBam(scratchBam, param=param);toc();
	message("\nCalculating density...\t", appendLF=FALSE);
	name <- names(sb)
	#readsGR<-GRanges(seqnames=as.vector(sb[[name]]$rname),IRanges(start=sb[[name]]$pos,width=1),strand=sb[[name]]$strand);
	#wholeChrom = GRangesForUCSCGenome(genomeBuild, chrom=as.character(runValue(seqnames(singleChromGR))))
	#wholeChromData <- import.bw(bamFile, BigWigSelection(wholeChrom), asRangedData=F);
	#cov = coverage(readsGR) #resize(readsGR,1)

	#using IRanges:
	reads = IRanges(start=sb[[name]]$pos,width=1)
	cov = coverage(reads);

	if (useRle) { #if using GenomicRanges, sub "cov" with "cov[[1]]" to access the innards. If using IRanges, use "cov".
		completeViews = successiveViews(cov, width=rep((windowSize*2+1),length(cov)/windowStep), gapwidth=-((windowSize*2+1)-windowStep))
		windowMulsBroad = viewMuls(completeViews, na.rm=FALSE, kernelWeights);
		#object.size(completeViewsR)/1000/1000
	} else {
		completeViews = successiveViews(as.vector(cov), width=rep((windowSize*2+1),length(cov)/windowStep), gapwidth=-((windowSize*2+1)-windowStep))
		windowMulsBroad = viewMuls(completeViews, na.rm=FALSE, kernelWeights);
		#object.size(completeViews)/1000/1000
	}
	windowCount = length(completeViews);
	rm(completeViews); gc();
	densityEstimate = signif(windowMulsBroad, digits=3); toc();
	rm(windowMulsBroad); gc();
	if (outputDensity || genomeWideCutoff) {
		#If the user requests to save the the density to disk, save here. Otherwise they will be discarded.
		#If the user requests a genome-wide cutoff, even without wanting to save the output, we 
		#store COMPLETE densities temporarily, until calculating the genome-wide cutoff, reading the density in again, and calling peaks.
		#Then we delete the temp files.
		message("\nStoring density...\t", appendLF=FALSE);
		wiggleFileName = paste(scratchDir, chrom, ".density.", windowSize, "-", windowStep, ".wig", sep="")
		write(paste("fixedStep  chrom=", chrom, " start=",windowSize+1,"  step=", windowStep, "  span=", 1, sep=""), file=wiggleFileName, sep="\t")
		#tic(); write.table(densityEstimates2, row.names=FALSE, col.names=FALSE, quote=FALSE, file=wiggleFileName, append=TRUE); toc();
		cwrite(as.vector(densityEstimate), as.numeric(length(densityEstimate)), file=wiggleFileName, showWarnings=FALSE); 

	}
	if (genomeWideCutoff) {
		#If the user requests a genome-wide cutoff, then we also will spit out a sample (default: 10%) 
		#of every chromosome into a temporary file. This is faster than using all the data, 
		#and more accurate than using only a single chromosome, which may be biased.
		sampleFileName = paste(scratchDir, chrom, ".sample.", windowSize, "-", windowStep, ".wig", sep="")
		file.create(sampleFileName);
		cwriteStep(as.vector(densityEstimate), as.numeric(length(densityEstimate)), file=sampleFileName, step=10, showWarnings=FALSE);
	}
	toc();
	if (!genomeWideCutoff) {
		#The alternative is to use a different cutoff for each chromosome. This can lead to problems if
		#some chromosomes have strange distributions or less data.
		looseQuantile = quantile(densityEstimate[densityEstimate > 0], probs=c(.95))
		windowCoords = seq(from=51, by=5, length.out=windowCount)
		peakCoords = defineBroadPeaks(chrom, densityEstimate, windowCoords, scratchDir, looseQuantile,windowStep);
		defineNarrowPeaks(chrom, densityEstimate, windowCoords, scratchDir, peakCoords);
	}
	return(windowCount);
}
