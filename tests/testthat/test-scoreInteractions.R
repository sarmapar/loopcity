test_that("scores are generated", {
  hicfile <- system.file("extdata", "GM12878_chr22.hic",
                         package = "loopcity")

  mergedLoops <-
      GM12878_10KbLoops |>
      mariner::assignToBins(binSize = 10e3) |>
      mergeAnchors(1)
  GenomeInfoDb::seqlevelsStyle(mergedLoops) <- "ENSEMBL"

  connections <- connectLoopAnchors(mergedLoops, 1e6)

  scores <- scoreInteractions(x = connections, hicFile = hicfile, norm = "NONE",
                    loopCalls = mergedLoops)


  assignCommunities(InteractionSet::interactions(scores))

})


test_that("loops with multiple widths are handled correctly", {
    ## generate loops of various widths
    regions5k <- GenomicRanges::GRanges(seqnames = "22",
                                        IRanges::IRanges(seq(5e3,95e3,5e3),
                                                         width = 5001))

    regions10k <- GenomicRanges::GRanges(seqnames = "22",
                                         IRanges::IRanges(seq(10e3,90e3,10e3),
                                                          width = 10001))

    loops <- c(InteractionSet::GInteractions(anchor1 = c(1,3),
                                             anchor2 = c(7,9),
                                             regions = regions5k),
               InteractionSet::GInteractions(anchor1 = c(1,3),
                                             anchor2 = c(7,9),
                                             regions = regions10k))
})

