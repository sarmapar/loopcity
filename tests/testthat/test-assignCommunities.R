### Setup inputs
hicfile <- system.file("extdata", "GM12878_chr22.hic",
                       package = "loopcity")

loops <- GM12878_10KbLoops
GenomeInfoDb::seqlevelsStyle(loops) <- "ENSEMBL"

connections <- connectLoopAnchors(loops, 1e6)

scores <- scoreInteractions(x = connections,
                            hicFile = hicfile,
                            norm = "NONE",
                            pseudo = 5)

loops_with_scores <- InteractionSet::interactions(scores)

test_that("assignCommunities generates GInteractions", {
    assignCommunities(loops_with_scores) |>
        expect_s4_class("GInteractions")
})

test_that("Various cluster types can be used",{
    ## fast greedy
    assignCommunities(loops_with_scores,
                      clusterType = "fast_greedy") |>
        expect_s4_class("GInteractions")

    ## walktrap
    assignCommunities(loops_with_scores,
                      clusterType = "walktrap") |>
        expect_s4_class("GInteractions")

    ## infomap
    assignCommunities(loops_with_scores,
                      clusterType = "infomap") |>
        expect_s4_class("GInteractions")

    ## label prop
    assignCommunities(loops_with_scores,
                      clusterType = "label_prop") |>
        expect_s4_class("GInteractions")

    ## edge betweenness
    assignCommunities(loops_with_scores,
                      clusterType = "edge_betweenness") |>
        expect_warning("highest modularity")

    ## invalid cluster type
    assignCommunities(loops_with_scores,
                      clusterType = "random_mess") |>
        expect_error("not a valid clustering algorithm")
})

test_that("loops are checked correctly", {
    ## loops must be a GInteractions object
    test_df <- data.frame(x = c(1,2,3), y = c(3,6,9))
    test_char <- "test"

    assignCommunities(test_char) |>
        expect_error("GInteractions")

    assignCommunities(test_df) |>
        expect_error("convert a data-frame-like object")
})

test_that("scores are checked correctly", {
    ## `scores` must be a column in `loops` or a numeric vector
    ## of the same length
    test_scores <- c(2,4,6,8,10)
    assignCommunities(loops_with_scores, score = test_scores) |>
        expect_error("must be the same length")

    test_scores_char <- "test"
    assignCommunities(loops_with_scores, score = test_scores_char) |>
        expect_error("must be a numeric vector")

    ## Providing scores overwrites `score` column in `loops`
    test_scores_correct <- 1:length(loops_with_scores)
    assignCommunities(loops_with_scores, score = test_scores_correct) |>
        expect_s4_class("GInteractions")

    ## Error if no scores provided or in `loops`
    loops_without_scores <- connections
    assignCommunities(loops_without_scores) |>
        expect_error("does not contain a column")

    ## Works properly is no scores in `loops` but provided separately
    ## set seed since clustering involves some randomness
    set.seed(555)
    scoresInLoopsOutput <- assignCommunities(loops_with_scores)
    scoresSeperateOutput <- assignCommunities(loops_without_scores,
                                              scores = as.numeric(
                                                  loops_with_scores$score))

    expect_identical(scoresInLoopsOutput$loopCommunity,
                     scoresSeperateOutput$loopCommunity)
})

test_that("leiden_resolution checks work properly", {
    ## Give a warning if leiden_resolution is provided but not used
    assignCommunities(loops_with_scores,
                      clusterType = "fast_greedy",
                      leiden_resolution = 0.5) |>
        expect_warning("will not be used")

    ## Give a warning if more than one resolution is provided
    assignCommunities(loops_with_scores,
                      leiden_resolution = c(1,2,3)) |>
        expect_warning("first will be used")

    ## Error if leiden_resolution is not numeric
    assignCommunities(loops_with_scores,
                      leiden_resolution = "test") |>
        expect_error("must be numeric")
})
