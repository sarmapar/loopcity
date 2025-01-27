plotCommunities <- function(loops, communityNums, buffer, hicFile){
    # params
    regions <- InteractionSet::regions(loops) |>
        dplyr::filter(anchorCommunity %in% communityNums)

    ### for future: one page per chrom, right now error
    windowChrom <- unique(GenomeInfoDb::seqnames(regions))
    if(length(chrom) > 1){
        abort("multi-chrom plotting not currently supported")
    }

    commStart = min(IRanges::start(regions))
    commEnd = max(IRanges::end(regions))

    windowStart = commStart - buffer
    windowEnd = commEnd

    # filter loops to window & scores greater than 0 (also removes NAs)
    windowLoops_i <- GenomicRanges::findOverlaps(regions, loops) |>
        S4Vectors::to()
    windowLoops <- loops[windowLoops_i]
    windowLoops <- windowLoops[which(as.logical(windowLoops$score > 0))]

    z = plotgardener::readHic(file = hicFile,
                chrom = windowChrom,
                chromstart = windowStart,
                chromend = windowEnd,
                resolution = 10e3)

    if(nrow(z) > 0){
        zmax = as.numeric(ceiling(quantile(z[z[,3] >= 10,3])[4])) + 200
    }

    # make plot
    plotgardener::pageCreate(width = 10, height = 7, showGuides = F)

    hicPlot <- plotgardener::plotHicRectangle(data = hicFile,
                                chrom = windowChrom,
                                chromstart = windowStart,
                                chromend = windowEnd,
                                resolution = 10e3,
                                zrange = c(0,100),
                                x = 0.5, y = 2.5, width = 9, height = 2.5)

    plotgardener::annoGenomeLabel(hicPlot, x = 0.5, y = 6.3, scale = "Mb")
    plotgardener::annoHeatmapLegend(hicPlot, x = 0.25, y = 2.5,
                                    width = 0.2, height = 1.8)

    anchors <- regions(windowLoops)
    plotRanges(data = anchors,
               chrom = windowChrom,
               collapse = T,
               chromstart = windowStart-buffer,
               chromend = windowEnd+buffer,
               x = 0.5, y = 5,
               height = 0.25, width = 9,
               fill = colorby("anchorCommunity",
                              palette = colorSkipper, scalePerRegion = TRUE))

    # find max score
    max_score <- max(windowLoops$score)

    #plot arches not in a community
    loopsNoCommunity <- windowLoops[which(windowLoops$loopCommunity==0)]
    if(length(loopsNoCommunity) > 0){
        arches <- plotPairsArches(loopsNoCommunity, chrom = windowChrom,
                                  chromstart = windowStart-buffer,
                                  chromend = windowEnd+buffer,
                                  flip = TRUE,
                                  archHeight = "score",
                                  range = c(0, max_score),
                                  fill = "grey",
                                  x = 0.5, y = 5, width = 9, height = 1.25)
    }

    #plot arches in a community over top
    loopsInCommunity <- windowLoops[which(windowLoops$loopCommunity>0 & windowLoops$score > 0)]
    if(length(loopsInCommunity) > 0){
        arches2 <- plotPairsArches(loopsInCommunity, chrom = windowChrom,
                                   chromstart = windowStart-buffer,
                                   chromend = windowEnd+buffer,
                                   flip = TRUE,
                                   archHeight = "score",
                                   range = c(0, max_score),
                                   fill = colorby("loopCommunity",
                                                  palette = viridis),
                                   alpha = 0.4,
                                   x = 0.5, y = 5, width = 9, height = 1.25)
    }
}
