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
