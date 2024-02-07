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
        ## decide how to choose best pseudocount, ROC only?
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
        roc <- pROC::roc(predictor = as.numeric(inters$score), response = inters$called_loop)
        auc <- pROC::auc(roc)
        pseudoAUCs <- dplyr::bind_rows(pseudoAUCs,
                                       data.table::data.table(
                                           auc = as.numeric(auc),
                                           psc = pseudo))
    }

    ## return pseudocounts with maximum AUC
    return(pseudoAUCs$psc[which.max(pseudoAUCs$auc)])
}
