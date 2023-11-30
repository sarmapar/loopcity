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
connectLoopAnchors <- function(loops, overlapDist){
    anchors <- InteractionSet::reduceRegions(loops) |>
        InteractionSet::regions()

    overlaps <- anchors %>%
        IRanges::findOverlaps(maxgap = overlapDist) %>%
        data.table::as.data.table() %>%
        dplyr::filter(queryHits != subjectHits)

    newInteractions <- InteractionSet::GInteractions(
        anchor1 = anchors[overlaps$queryHits],
        anchor2 = anchors[overlaps$subjectHits]) |>
        InteractionSet::swapAnchors(mode = "order") |>
        unique()

    return(newInteractions)
}
