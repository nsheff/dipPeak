
#This function just collects 2 functions into one:
#First, defines Broad peaks,
#Second, defines narrow peaks.
defineAllPeaks = function(chrom, scratchDir, windowSize, windowStep, windowCount, looseQuantile) {
		cat("\nChrom:", chrom, "\t");
		wiggleFileName = paste(scratchDir, chrom, ".density.", windowSize, "-", windowStep, ".wig", sep="")
		densityEstimate = scan(wiggleFileName, skip=1, allowEscapes=TRUE)
		windowCoords = seq(from=(windowSize+1), by=windowStep, length.out=windowCount[[chrom]])
		peakCoords = defineBroadPeaks(chrom, densityEstimate, windowCoords, scratchDir, looseQuantile, windowStep);
		if(!is.null(peakCoords))
		defineNarrowPeaks(chrom, densityEstimate, windowCoords, scratchDir, peakCoords);
	}

