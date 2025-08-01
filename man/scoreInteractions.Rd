% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scoreInteractions.R
\name{scoreInteractions,GInteractions,character-method}
\alias{scoreInteractions,GInteractions,character-method}
\alias{scoreInteractions,InteractionArray,missing-method}
\title{Score loop counts using enrichment of foreground over background}
\usage{
\S4method{scoreInteractions}{GInteractions,character}(
  x,
  hicFile,
  fgSize = 0,
  bg = "bowtie",
  bgGap = 0,
  bgSize = 5,
  pseudo = 5,
  pruneUnder = 0,
  truncate = F,
  resolution,
  norm = "NONE"
)

\S4method{scoreInteractions}{InteractionArray,missing}(x, fgSize, bg, bgGap, bgSize, pseudo, pruneUnder, truncate)
}
\arguments{
\item{x}{GInteractions object or InteractionArray object containing
interactions of interest
 See also `mariner::pixelsToMatrices`}

\item{fgSize}{manhattan distance surrounding center pixel to
include in foreground selection
See also parameter `n` in `mariner::selectRadius`}

\item{bg}{character keyword for background selection shape
One of `c("bowtie", "donut")`}

\item{bgGap}{number of pixels between foreground and background shapes}

\item{bgSize}{Integer describing the width (in pixels) of the background
shape}

\item{pseudo}{number of pseudo counts to add to raw Hi-C counts}

\item{pruneUnder}{numeric, all loops with raw counts less than this value
will have score set to NA}

\item{truncate}{logical indicating whether to set any values above
the 99th percentile to the 99th percentile value}

\item{resolution}{Integer (numeric) describing the resolution (range widths)
 of the paired data. Used only if GInteractions object is supplied for x.
By default, will match the width of regions in `x`}

\item{norm}{String (length one character vector)
describing the Hi-C normalization to apply. Use
`strawr::readHicNormTypes()` to see accepted values
for each file in `files`.
Required only if GInteractions object is supplied for x.}

\item{file}{Character file paths to a `.hic` file. Required only if
GInteractions object is supplied for x.
See also `qc_scoreInteractions`}

\item{interactions}{GInteraction object of interactions from Hi-C data}
}
\description{
Given a hi-c file or InteractionArray, pseudocounts are added to raw counts
and the median of the foreground over the background is added as a score for
each given interaction.
}
\examples{
hicFile <- "inst/extdata/GM12878_chr22.hic"

mergedLoops <- mergeAnchors(GM12878_10KbLoops, 1)
connections <- connectLoopAnchors(mergedLoops, 1e6)

scoreInteractions(connections, hicFile, mergedLoops)

}
