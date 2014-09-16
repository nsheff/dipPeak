#CWriteR
install.packages("~/Rncs/CWriteR", repos=NULL)
library(CWriteR)
a = rnorm(1000)
cwrite(a, length(a), "output.txt")


#IRangeKernels
install.packages("~/Rncs/IRangeKernels", repos=NULL)
library(IRangeKernels)
a = successiveViews(Rle(as.integer(rep(5, 50))), c(10, 10, 10))
class(a)
y = 1:10 + .001
viewSums(a)

viewMuls(a, y=as.numeric(y))
a

#dipPeak
install.packages("~/Rncs/dipPeak", repos=NULL)
library(dipPeak)
dipPeaks(bamFile="~/fhgfs/share/data/encodeRegions.bam", bigWigOut="encode.bw", outDir="hello", limitChrom=c("chr22", "chr5"), chromInfo="~/fhgfs/share/data/ucsc/chromInfo.hg19.txt")

