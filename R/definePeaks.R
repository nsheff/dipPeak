
#This function just collects 2 functions into one:
#First, defines Broad peaks,
#Second, defines narrow peaks.
defineAllPeaks = function(chrom, SCRATCHDIR, windowSize, windowStep, windowCount, looseQuantile) {
		cat("\nChrom:", chrom, "\t");
		wiggleFileName = paste(SCRATCHDIR, chrom, ".density.", windowSize, "-", windowStep, ".wig", sep="")
		densityEstimate = scan(wiggleFileName, skip=1, allowEscapes=TRUE)
		windowCoords = seq(from=51, by=5, length.out=windowCount[[chrom]])
		peakCoords = defineBroadPeaks(chrom, densityEstimate, windowCoords, SCRATCHDIR, looseQuantile, windowStep);
		if(!is.null(peakCoords))
		defineNarrowPeaks(chrom, densityEstimate, windowCoords, SCRATCHDIR, peakCoords);
	}

