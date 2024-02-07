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
            bowtie = mariner::selectTopLeft(n = n, buffer = buffer) +
                mariner::selectBottomRight(n = n, buffer = buffer),
            donut = mariner::selectOuter(n = n, buffer = buffer)
    )
}

## Generic for scoreInteractions
#' @rdname calcLoopEnrichment
#' @export
setGeneric("scoreInteractions",
           function(x,
                    hicFile,
                    loopCalls,
                    fgSize = 0,
                    bg = "bowtie",
                    bgGap = 0,
                    bgSize = 5,
                    pseudo,
                    pruneUnder = 0,
                    truncate = F,
                    resolution,
                    norm = "KR")
               standardGeneric("scoreInteractions"))


#' Internal for scoreInteractionsFromMats
#' @inheritParams scoreInteractions
#' @importFrom mariner calcLoopEnrichment
#' @noRd
#' @keywords internal
.scoreInteractionsFromMats <- function(x,loopCalls,fgSize,bg,bgGap,bgSize,
                                       pseudo,pruneUnder,truncate){

    ## Get buffer from counts matrices
    buffer <- (dim(BiocGenerics::counts(x))[1] - 1)/2

    ## Check that buffer is large enough for bg and fg
    if(buffer < fgSize + bgGap + bgSize){
        abort(c("`x` dimensions are too small",
                "x"=glue('`dim(counts(x))` must be at least {fgSize+bgGap+bgSize}'),
                "i"=glue('`dim(counts(x))` is {dim(counts(x))}'),
                "i"="Try using a smaller value for fgSize, bgGap, or bgSize,
                or re-run with loop and hi-c files."))
    }

    ## Select center pixel/foreground
    cp_sel <- mariner::selectCenterPixel(mhDist = 1:fgSize, buffer = buffer)

    ## Select background based on bg keyword
    bg_sel <- .pickSelection(bg, n = bgSize, buffer = buffer)

    ## Calculate best pseudocounts if none are provided
    # allow for range of pseudocounts to do optimization
    if(missing(pseudo)){
        ## by default, use 10th-90th deciles of counts for potential
        ## pseudocounts
        pseudo <- BiocGenerics::counts(x) |>
            as.numeric() |>
            quantile(seq(.1,.9,.1), na.rm=T) |>
            unique()
    }

    if(length(pseudo) != 1){
        if(missing(loopCalls)){
            abort(c("`loopCalls` must be provided to optimize pseudocounts.",
                    "i"="Either provide a single value for `pseudo`
                    or provide a GInteractions object for `loopCalls`."))
        }

        if(length(pseudo) > 1){
            pseudoValues <- pseudo
            pseudo <- .optimizePseudocounts(values = pseudoValues,
                                            mats = x,
                                            fgSize = fgSize,
                                            bg = bg,
                                            bgGap = bgGap,
                                            bgSize = bgSize,
                                            loopCalls = loopCalls,
                                            pruneUnder = pruneUnder,
                                            truncate = truncate)
        }

        message(glue("{pseudo} pseudocounts added to raw counts,
                for quality control info use `qc_scoreInteractions()`"))

    }

    x$pseudocounts <- pseudo

    ## Function to divide median of fg and bg with pseudocounts
    scoringFunction <- function(fg,bg,pseudo=pseudo){
        median((fg + 1 + pseudo)/(bg + 1 + pseudo), na.rm = T)
    }

    ## Calculate enrichment scores
    scores <- mariner::calcLoopEnrichment(x = x, fg = cp_sel, bg = bg_sel,
                                          FUN = scoringFunction)

    ## prune low scores
    scores[which(as.logical(scores<=pruneUnder)),] <- 0

    ## truncate high scores to 99th percentile
    if(truncate){
        newMax <- quantile(scores,.99)
        scores[which(scores>newMax),] <- newMax
    }

    GenomicRanges::mcols(x)$score <- scores

    return(x)
}

#' Internal for scoreInteractionsFromMats
#' @inheritParams scoreInteractions
#' @importFrom mariner calcLoopEnrichment
#' @noRd
#' @keywords internal
.scoreInteractionsFromFiles <- function(x,hicFile,loopCalls,fgSize,bg,bgGap,
                                        bgSize,pseudo,pruneUnder,truncate,
                                        resolution,norm){

    ## get resolution from loops
    if(missing(resolution)){
        uniqueWidths <- InteractionSet::regions(x) |>
            IRanges::width() |>
            unique()
        if(length(uniqueWidths) > 1){
            abort(c("Regions in `x` must be the same resolution.",
                    "x"=glue('`width(regions(x))` include
                             {length(uniqueWidths)} different widths.'),
                    "i"=glue("Providing a value for `resolution` will put
                             regions into equal sized bins.")))
        }
        resolution = uniqueWidths-1
    }

    ## pull hi-c contacts
    ## Set params
    buffer = fgSize + bgGap + bgSize

    ## Extract matrices around each loop
    mats <- mariner::pixelsToMatrices(x=x, buffer=buffer) |>
        mariner::pullHicMatrices(
            files=hicFile,
            binSize=resolution,
            half="upper",
            norm = norm
        )

    ## Remove buffer from anchors
    InteractionSet::interactions(mats) <- x

    ## Pass in parameters, including optional ones, to the next function
    args <- list(x = mats,
                 fgSize = fgSize,
                 bg = bg,
                 bgGap = bgGap,
                 bgSize = bgSize,
                 pruneUnder = pruneUnder,
                 truncate = truncate)
    if (!missing(loopCalls)){
        args$loopCalls <- loopCalls
    }
    if (!missing(pseudo)){
        args$pseudo <- pseudo
    }

    do.call(".scoreInteractionsFromMats", args)
}

#' Score loop counts using enrichment of foreground over background
#'
#' Given a hi-c file or InteractionArray, pseudocounts are added to raw counts
#' and the median of the foreground over the background is added as a score for
#' each given interaction.
#'
#' @param x GInteractions object or InteractionArray object containing
#' interactions of interest
#'  See also `mariner::pixelsToMatrices`
#' @param file Character file paths to a `.hic` file. Required only if
#'  GInteractions object is supplied for x.
#' @param loopCalls GInteractions object containing "truth" loop calls for
#' optimizing pseudocounts
#'  See also `qc_scoreInteractions`
#' @param fgSize manhattan distance surrounding center pixel to
#'  include in foreground selection
#'  See also parameter `n` in `mariner::selectRadius`
#' @param bg character keyword for background selection shape
#'  One of `c("bowtie", "donut")`
#' @param bgGap number of pixels between foreground and background shapes
#' @param bgSize Integer describing the width (in pixels) of the background
#' shape
#' @param pseudo Integer or integer vector of pseudocounts to add to raw counts.
#'  If more than one value is provided, the optimal value will be determined
#' and used.
#'  If no value is provided, the optimal value between 1-200 will be determined
#' and used.
#'  To optimize pseudocounts, `loopCalls` must be provided.
#'  See also `qc_ScoreInteractions`
#' @param pruneUnder numeric, all scores under this value will be set to 0
#' @param interactions GInteraction object of interactions from Hi-C data
#' @param pseudo number of pseudo counts to add to raw Hi-C counts
#' @param truncate logical indicating whether to set any values above
#'  the 99th percentile to the 99th percentile value
#' @param resolution Integer (numeric) describing the resolution (range widths)
#'  of the paired data. Used only if GInteractions object is supplied for x.
#' By default, will match the width of regions in `x`
#' @param norm String (length one character vector)
#'  describing the Hi-C normalization to apply. Use
#'  `strawr::readHicNormTypes()` to see accepted values
#'  for each file in `files`.
#'  Required only if GInteractions object is supplied for x.
#'
#' @importFrom mariner calcLoopEnrichment
#' @rdname scoreInteractions
#' @export
#'
#' @examples
#' ## add examples here
#'
setMethod("scoreInteractions",
          signature(x="GInteractions",
                    hicFile="character"),
          definition=.scoreInteractionsFromFiles)


#' Calculate loop enrichment over background.
#'
#' @rdname calcLoopEnrichment
#' @export
setMethod("scoreInteractions",
          signature(x="InteractionArray",
                    hicFile="missing"),
          definition=.scoreInteractionsFromMats)

