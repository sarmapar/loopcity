regions5k <- GenomicRanges::GRanges(seqnames = "chr1",
                                    IRanges::IRanges(seq(5e3,95e3,5e3),
                                                     width = 5001))

test_that("overlaps the correct distance", {
    inters <- InteractionSet::GInteractions(anchor1 = c(1,3,5),
                                            anchor2 = c(8,10,8),
                                            regions = regions5k)

    connectLoopAnchors(inters, 11e3)

    expected <- InteractionSet::GInteractions(anchor1 = c(1,3,5,8),
                                              anchor2 = c(3,5,8,10),
                                              regions = regions5k) |>
        InteractionSet::reduceRegions()

    connectLoopAnchors(inters, 11e3) |>
        expect_identical(expected)
})
