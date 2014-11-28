	# provide extremaS matrix, and x = c(x1,x2), two areas to search between,
	# will return the lowest minima between those two extrema.
findDivs = function(ind, extremaS) {
	x1=ind[1]
	x2=ind[2]
	res = which(extremaS[x1:x2,3] == min(extremaS[x1:x2,3]))
	return(x1+res[1]-1);
}


#Given some data and a quantile, fit a gamma distribution to the data
#And return the value of quantile for that paramaterized distribution.
fitGammaCutoff = function(d, quant=0.9) {
	v = var (d);
	m = mean(d);
	theta = v/m;
	kappa = m/theta;
	cutoff = qgamma(quant, shape=kappa, scale=theta)
	return(cutoff);
}

	#####################################################
	#Define Broad peaks
	#####################################################
defineBroadPeaks = function(chrom, densityEstimate, windowCoords, scratchDir, looseQuantile,windowStep) {
	cat("Defining broad peaks...\t");
	peakCalls = Rle(densityEstimate>looseQuantile)
	#restrict peak calls to only those that span a certain distance without dropping below the threshold.
	rTrueCumLen = cumsum(runLength(peakCalls))
	choose = (which(runValue(peakCalls) & runLength(peakCalls) > 100/windowStep)-1)
	if (length(choose) <1 ) { #found no broad peaks 
		return(NULL);
	}
	#pick out the starts and ends of the peaks we've selected
	peakStarts = rTrueCumLen[choose]
	peakEnds = rTrueCumLen[choose]+runLength(peakCalls)[choose+1]
	length(peakStarts) #how many peaks ?
	options(scipen=15) #don't print scientific notation
	highScores = sapply(mapply(seq, peakStarts, peakEnds), function(x) { return(max(densityEstimate[x])); } )
	#convert scores to the UCSC bed format 1-1000 scale.
	highScores = round( (highScores/max(highScores)) * 1000 ); #signif(highScores,4))
	write(file=paste(scratchDir, "/", chrom, ".peaks.b.bed", sep=""), rbind(chrom, windowCoords[peakStarts], windowCoords[peakEnds], 1:length(peakStarts), highScores, ncol=5, sep="\t"); toc();
	peakCoords = cbind(peakStarts, peakEnds)
	return(peakCoords);
}
	#####################################################
	#Define Narrow peaks
	#####################################################
defineNarrowPeaks = function(chrom, densityEstimate, windowCoords, scratchDir, peakCoords) {
	cat("Defining narrow peaks...\t",chrom);
	peakList = apply(peakCoords, 1, findPeaks, densityEstimate, windowCoords)
	peaks = do.call(rbind, peakList)
	#convert scores to the UCSC bed format 1-1000 scale.
	peaks[,3] = round( ( peaks[,3]/max( peaks[,3])) * 1000 );
	write(file=paste(scratchDir, "/", chrom, ".peaks.dips.bed", sep=""), rbind(chrom, windowCoords[peaks[,1]], windowCoords[peaks[,2]], 1:nrow(peaks), peaks[,3]),  ncol=5, sep="\t"); toc();
}




pkFig = function(extremaS) {
	val = densityEstimates2[peakStart:peakEnd]
	par(mfrow=c(3,1), mar=c(2,3,2,1))
	plot(val, type="l", lwd=2, main="signal")
	abline(v=deriv1Zeros, col=rgb(.6,.6,.9))
	abline(v=tips, col="magenta");
	points(extremaS[extremaS[,2]==-1,1], extremaS[extremaS[,2]==-1,3], bg="red", pch=25, cex=2)
	points(extremaS[extremaS[,2]==1,1], extremaS[extremaS[,2]==1,3], bg="darkgreen", pch=24, cex=2)
	points(divMatrix, rep(max(peakValues),length(divMatrix)), lwd=2)
	plot(deriv1, type="l", col="blue", main="1st derivative")
	abline(h=0, col=rgb(.8,.8,.8))
	abline(v=deriv1Zeros, col=rgb(.6,.6,.9))
	plot(deriv2, type="l",  col="red", main="2nd derivative")
	abline(h=0, col=rgb(.8,.8,.8))
	abline(v=deriv1Zeros, col=rgb(.6,.6,.9))
	#abline(v=deriv1Zeros[deriv2[deriv1Zeros] <=0], col="cyan")
}


killdip = function(zen, extremaS) {
	if (nrow(extremaS) <4) { return(NULL); }
	if (zen == 2) { return(3); }
	if (zen == nrow(extremaS)-1 ) { return( nrow(extremaS)-2); }
	#if (extremaS[zen+2,3] > extremaS[zen-2,3]) { return (zen+1); } #merge this peak with the larger neighboring peak
	if (abs(extremaS[zen+2,1]-extremaS[zen,1]) < abs(extremaS[zen-2,1]-extremaS[zen,1])) { return (zen+1); }
	return(zen-1);
}

symdiff = function(x) {
	d = diff(c(x, rep(x[length(x)], 1)) , lag=1) + 
	diff(c(x, rep(x[length(x)], 2)) , lag=2) +
	diff(c(x, rep(x[length(x)], 3)) , lag=3) +
	diff(c(rep(x[1], 1), x) , lag=1) + 
	diff(c(rep(x[1], 2), x) , lag=2) +
	diff(c(rep(x[1], 3), x) , lag=3);
	return(d/6);
}

#Function takes a vector with coordinates of the start and end of a region,
#and returns a vector of the tips in that region
#depending on some tweakable parameters.
showfig=TRUE;
findPeaks = function(broadPeakCoords, densityEstimate, windowCoords, showfig=FALSE) {
	#cat(broadPeakCoords[1], "\t"); 	#debug information
	peakStart = broadPeakCoords[1];
	peakEnd = broadPeakCoords[2];
	peakWidth = windowCoords[peakEnd] - windowCoords[peakStart]
	peakValues = densityEstimate[peakStart:peakEnd]
	mpv = mean(peakValues)
	deriv1 = symdiff(peakValues) 	#first derivative of signal
	deriv2 = symdiff(deriv1)		#second derivative of signal
	signalRle = Rle(peakValues >= 0)								
	deriv1Rle = Rle(deriv1 >= 0)
	deriv1Zeros = cumsum(runLength(deriv1Rle))

	tips = deriv1Zeros[deriv2[deriv1Zeros] < -mpv*.01]; #peakValues[deriv1Zeros] > strictQuantile & 
	dips = c(1,na.omit(deriv1Zeros[deriv2[deriv1Zeros] > mpv*.01]), length(peakValues))

	if (length(tips) == 0) { return(NULL); }
	#build a matrix with all the information about the tips and dips:
	extrema = rbind(cbind(na.omit(tips), 1, peakValues[tips], deriv1[tips], deriv2[tips]), cbind(na.omit(dips), -1, peakValues[dips], deriv1[dips], deriv2[dips]))
	extremaS = extrema[order(extrema[,1]),]

	#In a legitimate matrix, the order will go: dip, tip, dip, tip, dip...etc. 
	#If there are two dips or tips in a row, we should delete one.
	badExtrema = which(runLength(Rle(extremaS[,2])) > 1)
	if (length(badExtrema) > 0) { 
		#two or more of the same extrema in a row: we should only keep one, but which one? The best one.
		cs = cumsum(runLength(Rle(extremaS[,2])))
		runs = mapply(seq,  c(0, cs[-length(cs)])+1, cs, SIMPLIFY=FALSE)
		todelete = lapply(runs, function(x) { #keep only the highest tip or the lowest dip.
			keep = which(extremaS[x,3] == abs(max(extremaS[x,3] * extremaS[x,2])))[1];
			return(setdiff(x, x[keep]));
		})
		#eliminate bad extrema
		extremaS = extremaS[-unlist(todelete), ]
		#reconstitute the dips and tips vectors with the new extrema.
		dips = extremaS[extremaS[,2] == -1,1] 
		tips = extremaS[extremaS[,2] == 1,1]
	}
	#Now, we want to only keep the tips that have sufficient peakStrength,
	#which measures not the total peak height, but how much it rises above its local dips.
	extremaDeriv = diff(extremaS[,3])
	peakRise = (abs(extremaDeriv [ which (extremaS[,2] == 1)-1]) + abs(extremaDeriv [ which (extremaS[,2] == 1)]))/2
	peakStrength = peakRise / peakValues[tips]

	#delete any tip with low peakStrength
	#each tip has two dips, we want to delete one of them. Which one? Perhaps the furthest away; 
	#this is what function killdip() does.
	delete = peakStrength < .5
	i = order(peakStrength)[1]
	while (peakStrength[i] < .5) {
		if (length(peakStrength) == 1) { return(NULL); }
		#cat(peakStrength[i]); #debug information
		dipToDelete = killdip(which(extremaS[,2]==1)[i],extremaS)
		extremaS = extremaS[-c(which(extremaS[,2]==1)[i], dipToDelete), ]
		extremaDeriv = diff(extremaS[,3])
		dips = extremaS[extremaS[,2] == -1,1] 
		tips = extremaS[extremaS[,2] == 1,1]
		peakRise = (abs(extremaDeriv [ which (extremaS[,2] == 1)-1]) + abs(extremaDeriv [ which (extremaS[,2] == 1)]))/2
		peakStrength = peakRise / peakValues[tips]
		i = order(peakStrength)[1]
	}
	
	#now take the final set of peaks and convert it back to correct coordinates.
	divMatrix = cbind(dips[-length(dips)], dips[-1])
	peaksFound = peakStart + divMatrix
	if(showfig) { pkFig(extremaS); }
	return(cbind(peaksFound, peakValues[tips]))
}

