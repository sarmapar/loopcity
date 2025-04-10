#' Helper function to optimize pseudocounts
#' @param x character vector of length 1 describing keyword
#'  (one of "bowtie", "donut")
#' @param n Integer describing the number of outer pixels to select.
#'  Must be length of one. See also `?mariner::selectRadius`
#' @returns a MatrixSelection object corresponding to the keyword
#'
#' @noRd
#' @keywords internal
.optimizePseudocounts <- function(mats,values,loopCalls,fgSize,bg,bgGap, bgSize,
                                  pseudo,pruneUnder,truncate){

    pseudoAUCs <- data.table::data.table(auc = numeric(), psc = numeric())
    for(pseudo in values){
        loopScores <- scoreInteractions(
            pseudo = pseudo,
            x = mats,
            fgSize = fgSize,
            bg = bg,
            bgGap = bgGap,
            bgSize = bgSize,
            pruneUnder = pruneUnder,
            truncate = truncate)

        ## pull out interactions from InteractionArray
        inters <- InteractionSet::interactions(loopScores)

        ## add column if interactions are present in original loop calls
        inters$called_loop <- FALSE
        inters$called_loop[inters %in% loopCalls] <- TRUE

        ## calculate area under the roc curve
        roc <- pROC::roc(predictor = as.numeric(inters$score),
                         response = inters$called_loop)
        auc <- pROC::auc(roc)
        pseudoAUCs <- dplyr::bind_rows(pseudoAUCs,
                                       data.table::data.table(
                                           auc = as.numeric(auc),
                                           psc = pseudo))
    }

    ## return pseudocounts with maximum AUC
    return(pseudoAUCs$psc[which.max(pseudoAUCs$auc)])
}



#' Helper function to determine if the node on the border of a community should
#' be assigned to both its current community and the neighboring community
#' @param loops GInteractions object containing all loops of interest, must also
#' contain a `score` metadata column
#' @param borderNode GRanges object of the border node to be analyzed, must also
#' contain a `anchorCommunity` metadata column
#' @param left logical indicating if borderNode is on the left border of its
#' assigned community
#' @returns GRanges object for borderNode with new community assignment
#'
#' @noRd
#' @keywords internal
.compareNodeToComm <- function(loops, borderNode, left){
    ## Suppress NSE notes in R CMD check
    queryHits <- loopCommunity <- comm <- score <- NULL

    ## Find loops that contain the borderNode (overlaps 1) and anchors in
    ## the neighboring community (overlaps 2)
    if(left){
        overlaps1 <- InteractionSet::findOverlaps(
            InteractionSet::anchors(loops,"second"), borderNode)
        overlaps2 <- InteractionSet::findOverlaps(
            InteractionSet::anchors(loops,"first"),
            InteractionSet::regions(loops)[which(
            InteractionSet::regions(loops)$anchorCommunity ==
                borderNode$anchorCommunity-1)])
    } else {
        overlaps1 <- InteractionSet::findOverlaps(
            InteractionSet::anchors(loops,"first"), borderNode)
        overlaps2 <- InteractionSet::findOverlaps(
            InteractionSet::anchors(loops,"second"),
            InteractionSet::regions(loops)[which(
            InteractionSet::regions(loops)$anchorCommunity ==
                borderNode$anchorCommunity+1)])
    }

    # Keep only interactions between border node and neighboring community nodes
    interactions_to_keep <- overlaps1[
        S4Vectors::queryHits(overlaps1) %in% S4Vectors::queryHits(overlaps2)]

    ## Compare the median score of all interactions within the neighboring
    ## community (comm_med) to all interactions between the border node and
    ## anchors in the neighboring community (node_med)
    filtered_scores <- loops$score[queryHits(interactions_to_keep)]
    node_med <- filtered_scores |> median()

    comm_med <- loops |>
        as.data.frame() |>
        dplyr::filter(loopCommunity == comm) |>
        dplyr::pull(score) |>
        median()

    ## If node_med >= comm_med, assign the border node to both communities
    if(!is.na(node_med) & node_med >= comm_med){
        if(left){
            borderNode$anchorCommunity <- list(c(comm-1,
                                                 borderNode$anchorCommunity))
        } else {
            borderNode$anchorCommunity <- list(c(borderNode$anchorCommunity,
                                                 comm+1))
        }
    }

    return(borderNode)
}
