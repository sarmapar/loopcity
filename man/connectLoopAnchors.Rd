% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/connectLoopAnchors.R
\name{connectLoopAnchors}
\alias{connectLoopAnchors}
\title{connectLoopAnchors}
\usage{
connectLoopAnchors(loops, overlapDist)
}
\arguments{
\item{loops}{GInteractions object of loops}

\item{overlapDist}{an integer representing the farthest distance in base pairs
to make a connection}
}
\value{
GInteractions object
}
\description{
connectLoopAnchors
}
\examples{
mergedLoops <- mergeAnchors(GM12878_10KbLoops, 1)

connectLoopAnchors(mergedLoops, 1e6)
}
