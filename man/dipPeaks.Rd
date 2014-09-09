\name{dipPeaks}
\alias{dipPeaks}
\title{Estimate density and find peaks}
\usage{
  dipPeaks(BAMFILE, INDEXFILE, SCRATCHDIRBASE, BIGWIGOUT,
    OUTDIR, CORES = 1, USE_RLE = TRUE,
    GENOME_WIDE_CUTOFF = TRUE, OUTPUT_DENSITY = TRUE,
    windowSize = 50, windowStep = 5)
}
\arguments{
  \item{BAMFILE}{Input file in bam format}

  \item{INDEXFILE}{Index file made from "samtools index" on
  the bam input file}

  \item{SCRATCHDIRBASE}{Directory for scratch output; use
  local disk for optimal processing. dipPeak writes a
  significant amount of data to disk as it runs, to reduce
  memory use. This I/O takes time, and will be faster if
  you use a local disk or fast mount spot.}

  \item{BIGWIGOUT}{Smoothed output density file in bigwig
  format (if you request an output density).}

  \item{OUTDIR}{Final folder (the files to keep will be
  moved here from the scratch dir).}

  \item{CORES}{Add additional cores (up to 23) for faster
  processing. (Default:1)}

  \item{USE_RLE}{Boolean; offers speed and memory
  improvement}

  \item{GENOME_WIDE_CUTOFF}{Boolean; If yes, calculates a
  genome-wide cutoff for peak detection. Otherwise, it does
  this on a per-chromosome basis}

  \item{OUTPUT_DENSITY}{Boolean; If yes, produces a bigwig
  density track.}

  \item{windowSize}{Choose the size of the smoothing window
  (default:50)}

  \item{windowStep}{Choose the window step (5 will take a
  window every 5 bases) (default:5)}
}
\description{
  Calculates a smoothed track from input count data, and
  identifies peaks in the data. Outputs the peaks and
  smoothed track to file.
}
\examples{
dipPeaks("sequences.bam", "sequences.bam.bai", "scratch", "density.bw", "OUTDIR", CORES=8);
dipPeaks("encodeRegions.bam", "encodeRegions.bam.bai", "scratch", "bw.bw", "finalOutput", 8)
}

