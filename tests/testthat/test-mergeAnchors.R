# create GRanges object to pull regions from
regions5k <- GenomicRanges::GRanges(seqnames = "chr1",
                   IRanges::IRanges(seq(5e3,95e3,5e3), width = 5001))

regions10k <- GenomicRanges::GRanges(seqnames = "chr1",
                      IRanges::IRanges(seq(10e3,90e3,10e3), width = 10001))


test_that("anchors merge to correct positions", {
    #create input ranges and interactions
    inters <- InteractionSet::GInteractions(anchor1 = c(2,1,3),
                                            anchor2 = c(5,7,7),
                                            regions = regions5k)

    #create expected output ranges and interactions with duplicates
    expectedDups <- InteractionSet::GInteractions(anchor1 = rep(2,3),
                                              anchor2 = c(5,7,7),
                                              regions = regions5k) |>
        InteractionSet::reduceRegions()

    #create expected output ranges and interactions without duplicates
    expectedNoDups <- InteractionSet::GInteractions(anchor1 = c(2,2),
                                                    anchor2 = c(5,7),
                                                    regions = regions5k) |>
        InteractionSet::reduceRegions()

    #generate outputs using mergeAnchors
    outputDups <- mergeAnchors(inters, pixelOverlap = 1, dropDups = F)
    mergeAnchors(inters, pixelOverlap = 1) |>
        expect_message("Duplicates dropped.")
    outputNoDups <- mergeAnchors(inters, pixelOverlap = 1)

    #compare expected output and function output
    expect_identical(outputDups, expectedDups)
    expect_identical(outputNoDups, expectedNoDups)
})

test_that("there is no message if there are no duplicates", {
    #create input interactions
    inters <- InteractionSet::GInteractions(anchor1 = c(2,1,3),
                                            anchor2 = c(5,7,9),
                                            regions = regions5k)

    #expect no message, dropDups is TRUE but there are no duplicates
    expect_no_message(mergeAnchors(inters, pixelOverlap = 1))

})

test_that("mixed range widths give an error when using pixelOverlap", {
    #create input interactions with mixed widths
    inters <- c(InteractionSet::GInteractions(anchor1 = c(1,3),
                                            anchor2 = c(7,9),
                                            regions = regions5k),
        InteractionSet::GInteractions(anchor1 = c(1,3),
                                      anchor2 = c(7,9),
                                      regions = regions10k))

    #expect an error about different widths
    expect_error(mergeAnchors(inters, pixelOverlap = 1),
                 "same width to use pixelOverlap")

})

