% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mergeAnchors.R
\name{mergeAnchors}
\alias{mergeAnchors}
\title{Merge neighboring loop anchor positions into one
 representative anchor.}
\usage{
mergeAnchors(loops, pixelOverlap, dropDups = T)
}
\arguments{
\item{loops}{GInteractions object of loops from Hi-C data}

\item{pixelOverlap}{numeric representing the number of
pixels away from an anchor is considered to be overlapping
For example, `pixelOverlap = 1` means only immediately neighboring pixels
will be merged and `pixelOverlap = 2` means all anchors within 2 pixels}

\item{dropDups}{logical determining whether duplicates should be dropped
from output}
}
\description{
Merge neighboring loop anchor positions into one
 representative anchor.
}
\examples{
mergeAnchors(loops = GM12878_10KbLoops, pixelOverlap = 1)

}
