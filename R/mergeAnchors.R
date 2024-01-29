#' Merge neighboring loop anchor positions into one
#'  representative anchor.
#' @param loops GInteractions object of loops from Hi-C data
#' @param pixelOverlap numeric representing the number of
#'  pixels away from an anchor is considered to be overlapping
#'  For example, `pixelOverlap = 1` means only immediately neighboring pixels
#'  will be merged and `pixelOverlap = 2` means all anchors within 2 pixels
#' @param dropDups logical determining whether duplicates should be dropped
#'  from output
#'
#' @export
#'
#' @examples
#' mergeAnchors(loops = GM12878_10KbLoops, pixelOverlap = 1)
#'

mergeAnchors <- function(loops, pixelOverlap, dropDups = T){
    ## Only include anchors used in the GInteractions object
    loops <- InteractionSet::reduceRegions(loops)

    ## Get resolution of loops
    resolution <-
        InteractionSet::regions(loops) |>
        IRanges::ranges() |>
        IRanges::width() |>
        unique()-1

    if(length(resolution) > 1){
        rlang::abort(c("Interaction ranges must have the same width ",
                       "to use pixelOverlap.",
                "x"=glue("There are {length(resolution)} different widths ",
                         "in `ranges(regions(loops))`."),
                "i"="Consider using bpOverlap instead."))
    }

    ## Create data table of all overlapping anchor indexes
    ## (one line per overlap)
    overlaps <- loops |>
        InteractionSet::regions() |>
        IRanges::findOverlaps(maxgap = pixelOverlap*resolution-2) |>
        data.table::as.data.table()

    ## Condense overlaps into list of all overlapping anchor indexes
    ## (one line per unique anchor)
    ## & give each chain a group number
    overlapGroups <- overlaps[, .(list(subjectHits)), by=queryHits] |>
        dplyr::rename(anchorId = queryHits, overlapIds = V1) |>
        dplyr::rowwise() |>
        dplyr::mutate(sameGroup = anchorId %in% dplyr::lead(overlapIds)) |>
        dplyr::ungroup() |>
        dplyr::mutate(group = cumsum(!sameGroup)) |>
        dplyr::select(-sameGroup)

    ## Add count column of how many times the anchor is in loop data
    overlapGroups$count <- InteractionSet::anchorIds(loops) |>
        unlist() |>
        as.numeric() |>
        sort() |>
        table() |>
        as.numeric()

    ## helper function to pick the middle value of a vector
    ## (or the largest of 2 middle vals)
    middle <- function(x) {x[length(x)/2 + 1]}

    ## For each group of anchors, collect all anchor IDs with the maximum counts (modes)
    ## If there is only one mode, that is the representative anchor
    ## If there is more than one mode, choose the middle-est one
    representativeAnchors <- overlapGroups |>
        dplyr::group_by(group) |>
        dplyr::mutate(repAnchorId = {modes <- anchorId[count == max(count)];
        ifelse(length(modes) > 1,
               middle(modes),
               modes)})

    ## Replace original loop anchor IDs with new anchor IDs and create
    ## new GInteractions object
    mergedLoops <- loops
    InteractionSet::regions(mergedLoops) <-
        InteractionSet::regions(loops)[representativeAnchors$repAnchorId]

    mergedLoops <- InteractionSet::GInteractions(
        anchor1 = InteractionSet::anchors(mergedLoops, 'first'),
        anchor2 = InteractionSet::anchors(mergedLoops, 'second'))

    ## Remove duplicates and give a message of how many interactions were
    ## dropped, if `dropDups` == TRUE
    if(dropDups){
        if(length(mergedLoops) != length(unique(mergedLoops))){
            message(glue("Duplicates dropped. {length(mergedLoops)} ",
                         "interactions reduced to {length(unique(mergedLoops))}",
                         " interactions."))
        }

        return(unique(mergedLoops))

    } else {
        return(mergedLoops)
    }
}
