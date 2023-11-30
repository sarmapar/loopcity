#' Helper function to return selection by keyword
#' @param x character vector of length 1 describing keyword
#'  (one of "bowtie", "donut")
#' @param n Integer describing the number of outer pixels to select.
#'  Must be length of one. See also `?mariner::selectRadius`
#' @returns a MatrixSelection object corresponding to the keyword
#' @noRd
#' @keywords internal
.pickSelection <- function(x = c("bowtie","donut"), n, buffer){
    x <- match.arg(x)
    switch (x,
            bowtie = selectTopLeft(n = n, buffer = buffer) +
                selectBottomRight(n = n, buffer = buffer),
            donut = selectOuter(n = n, buffer = buffer)
    )
}

#' Score loop count using enrichment of foreground over background
#' @param mats matrix of HiC count interactions.
#'  See also `mariner::pixelsToMatrices`
#' @param fgDist manhattan distance surrounding center pixel to
#'  include in foreground selection
#'  See also parameter `n` in `mariner::selectRadius`
#' @param bg character keyword for background selection
#'  One of `c("bowtie", "donut")`
#' @param pruneUnder numeric, all scores under this value will be set to 0
#' @param interactions GInteraction object of interactions from Hi-C data
#' @param pseudo number of pseudo counts to add to raw Hi-C counts
#' @param truncate logical indicating whether to set any values above
#'  the 99th percentile to the 99th percentile value
#'
#' @importFrom mariner calcLoopEnrichment
#' @rdname scoreInteractions
#' @export
#'
#' @examples
scoreInteractions<- function(mats,
                            fgDist = 0,
                            bg = "bowtie",
                            bgDist = 5,
                            pseudo = 0,
                            pruneUnder = 0,
                            truncate = F){

    ## Get buffer from counts matrices
    buffer <- (dim(counts(mats))[1] - 1)/2

    ## Select center pixel/foreground
    cp_sel <- selectCenterPixel(mhDist = 1:fgDist, buffer = buffer)

    ## Select background based on bg keyword
    bg_sel <- .pickSelection(bg, n = bgDist, buffer = buffer)

    ## Function to divide median of fg and bg with pseudocounts
    scoringFunction <- function(fg,bg,pseudo=pseudo){
        median((fg + 1 + pseudo)/(bg + 1 + pseudo), na.rm = T)
    }

    ## Calculate enrichment scores
    scores <- calcLoopEnrichment(x = mats, fg = cp_sel, bg = bg_sel, FUN = scoringFunction)

    ## prune low scores
    scores[which(scores<=pruneUnder),] <- 0

    ## truncate high scores
    if(truncate){
        newMax <- quantile(scores,.99)
        scores[which(scores>newMax),] <- newMax
    }

    return(scores)
}
