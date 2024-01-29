#' connectLoopAnchors
#'
#' @param loops GInteractions object of loops
#' @param overlapDist an integer representing the farthest distance in base pairs
#' to make a connection
#'
#' @return GInteractions object
#' @export
#'
#' @examples
#' mergedLoops <- mergeAnchors(GM12878_10KbLoops, 1)
#'
#' connectLoopAnchors(mergedLoops, 1e6)
connectLoopAnchors <- function(loops, overlapDist){

    ## Only include anchors used in the GInteractions object
    anchors <- InteractionSet::reduceRegions(loops) |>
        InteractionSet::regions()

    ## Find all overlaps within `overlapDist` base pairs,
    ## and remove anchors interacting with themselves
    overlaps <- anchors |>
        IRanges::findOverlaps(maxgap = overlapDist) |>
        data.table::as.data.table() |>
        dplyr::filter(queryHits != subjectHits)

    ## Create GInteractions object of connections
    ## Re-order so start2 is always larger than start1
    ## Remove duplicates
    newInteractions <- InteractionSet::GInteractions(
        anchor1 = anchors[overlaps$queryHits],
        anchor2 = anchors[overlaps$subjectHits]) |>
        InteractionSet::swapAnchors(mode = "order") |>
        unique()

    return(newInteractions)
}
