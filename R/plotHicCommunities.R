#' Plot Hi-C data with loops colored by community assignment.
#'
#' If a vector of starts and ends is provided, one page will be generated
#' per region.
#'
#' @param pdfName a character string giving the file path. See the `file`
#' parameter in `pdf()` for more information.
#' @param communities GInteractions object containing loops to plot and a
#' column called `loopCommunity` containing community assignments
#' @param hicFile path to .hic file to plot
#' @param chroms chromosome(s) of the region(s) to plot
#' @param starts starting genomic location(s) of the region(s) to plot
#' @param ends ending genomic location(s) of the region(s) to plot
#' @param zmax a numeric value specifying the maximum range of hic scores to
#' plot, where values greater than `zmax` will be set to`zmax`
#' @param colorPalette a function corresponding to a color palette
#'
#' @export
#'
#' @examples
#' add example
plotHicCommunities <- function(pdfName, communities, hicFile, norm = "SCALE",
                            chroms, starts, ends, zmax, colorPalette){
#### parameter checking --------------------------------------------------------
    ## if only one chromosome, replicate for each region
    if(length(chroms) == 1){
        chroms = rep(chroms, length(starts))
    }
    if(length(chroms) != length(starts) | length(chroms) != length(ends)){
        abort("`chroms`, `starts`, and `ends` must be the same length")
    }

    #todo check if colorPalette is a color palette
    if(missing(colorPalette)){
        colorPalette = loopcityColors
    }

#### plan plotgardener page ----------------------------------------------------
    plotW = 9
    hicH = 2
    labelH = 0.5
    archH = 0.75
    archSpace = 0.25
    border = 0.5

    pageHeight = border + hicH +
        archH + archSpace +
        labelH + border

    pageWidth = plotW + 2*border

    pdf(pdfName, width = pageWidth, height = pageHeight)

    for(i in 1:length(chroms)){
        # set region to current chr, start, end
        chr = chroms[i]
        start = starts[i]
        end = ends[i]

        # params
        regions <- InteractionSet::regions(communities)
        regions <-regions[
            GenomicRanges::seqnames(regions) == chr &
                GenomicRanges::start(regions) >= start &
                GenomicRanges::end(regions) <= end]

        if(missing(zmax)){
            z = plotgardener::readHic(file = hicFile,
                                      chrom = chr,
                                      chromstart = start,
                                      chromend = end,
                                      norm = norm,
                                      resolution = 10e3)

            if(nrow(z) > 0){
                zmax = as.numeric(ceiling(quantile(z[z[,3] >= 10,3])[4])) + 100
            } else{
                zmax = 0
                rlang::abort("`zmax` could not be calculated, please provide a value.")
            }
        }

        # make plot
        plotgardener::pageCreate(width = pageWidth,
                                 height = pageHeight,
                                 showGuides = F)
        # hic plot
        hicPlot <- plotgardener::plotHicRectangle(data = hicFile,
                                                  chrom = chr,
                                                  chromstart = start,
                                                  chromend = end,
                                                  norm = norm,
                                                  resolution = 10e3,
                                                  zrange = c(0,zmax),
                                                  x = border,
                                                  y = border,
                                                  width = plotW,
                                                  height = hicH)

        # filter for
        # - loops in region
        # - scores greater than 0 (also removes NAs)
        # - loops in a community
        windowLoops_i <- GenomicRanges::findOverlaps(regions,
                                                     communities) |>
            S4Vectors::to()
        windowLoops <- communities[windowLoops_i]
        windowLoops <- windowLoops[which(as.logical(windowLoops$score> 0) &
                                        lengths(windowLoops$loopCommunity)>0)]


        # if `source` column exists, annotate pixels
        # gray box = added loop, black box = original loop
        noSource = F
        if(is.null(windowLoops$source)){
            rlang::warn(glue::glue("Cannot differentiate added vs original ",
                                   "loops in annotated pixels or arches,",
                                   " no `source` column in `communities`."))

            tempLoops <- windowLoops
            S4Vectors::mcols(tempLoops) <- NULL
            plotgardener::annoPixels(hicPlot, tempLoops)

            noSource = T
        } else {
            # annotating original loops
            originalLoops <- windowLoops[which(as.logical(
                windowLoops$source == "original"))]
            S4Vectors::mcols(originalLoops) <- NULL

            plotgardener::annoPixels(hicPlot, originalLoops)

            # annotating added loops
            addedLoops <- windowLoops[which(as.logical(
                windowLoops$source == "added"))]
            S4Vectors::mcols(addedLoops) <- NULL

            if(length(addedLoops) > 0){
                plotgardener::annoPixels(hicPlot, addedLoops, col = "gray")
            }
        }

        # plot genome label and legend for hicplot
        plotgardener::annoGenomeLabel(hicPlot,
                                      x = border,
                                      y = pageHeight - border - labelH,
                                      scale = "Mb")

        plotgardener::annoHeatmapLegend(hicPlot, x = border-0.25, y = border,
                                        width = 0.2, height = hicH*(3/5))

        ## Plot arches, color by community and height based on log2(score)
        # Map unique values to colors
        currColors <- colorPalette(length(unique(windowLoops$loopCommunity)))
        value_to_color <- setNames(currColors,
                                   unique(windowLoops$loopCommunity))

        # Assign colors based on values
        windowLoops$color <- value_to_color[
            as.character(windowLoops$loopCommunity)]

        #calculate log 2 of scores for plotting arches
        windowLoops$score_lg <- log2(windowLoops$score)

        # find max score
        max_score <- max(windowLoops$score_lg)

        # plot arches
        if(max_score > 0){
        arches <- plotgardener::plotPairsArches(windowLoops,
                                                chrom = chr,
                                                chromstart = start,
                                                chromend = end,
                                                flip = TRUE,
                                                clip = TRUE,
                                                archHeight = "score_lg",
                                                range = c(0, max_score),
                                                x = border,
                                                y = border + hicH,
                                                height = archH,
                                                width = plotW,
                                                fill = windowLoops$color)
        }
    }
    dev.off()

    message(glue::glue("Plots created in "), pdfName)
}

