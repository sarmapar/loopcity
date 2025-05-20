#' Assign communities to interactions
#'
#' @param loops GInteractions object of interactions to include in the network
#' @param scores a numeric vector of scores to use as edge weights,
#' if `loops` contains a column named 'score', this column will be used
#' @param clusterType character vector of length 1 containing a keyword
#' for the clustering algorithm to use. Must be one of "fast_greedy",
#' "walktrap", "leiden", "infomap", "label_prop", or "edge_betweenness"
#' @param leidenResolution numeric between 0 and 1, higher resolutions lead to
#' more smaller communities, while lower resolutions lead to fewer larger
#' communities. See `resolution_parameter` in `igraph::cluster_leiden`
#' @param pruneUnder numeric, added loops with scores less than this value will
#' be removed
#'
#' @return GInteractions object with community metadata
#'
#' @importFrom rlang abort
#' @importFrom glue glue glue_collapse
#' @importFrom methods is
#' @importFrom stats density na.omit
#' @export
#'
#' @examples
#' hicFile <- "inst/extdata/GM12878_chr22.hic"
#'
#' mergedLoops <- mergeAnchors(GM12878_10KbLoops, 1)
#' connections <- connectLoopAnchors(mergedLoops, 1e6)
#' scores <- scoreInteractions(connections, hicFile, mergedLoops)
#'
#' assignCommunities(interactions(scores))
assignCommunities <- function(loops,
                              scores,
                              clusterType = "leiden",
                              leidenResolution = 0.1,
                              pruneUnder) {

    ## Suppress NSE notes in R CMD check
    anchorCommunity <- chr <- NULL

## Parameter checking ----------------------------------------------------------
    call_args <- match.call()

    ## Check if loops is a GInteractions object, with helpful message for
    ## converting data frame to GInteractions
    if(!is(loops,"GInteractions")){
        if(is(loops,"data.frame") | is(loops,"data.table") |
                                       is(loops,"DataFrame")){
            rlang::abort(c("`loops` must be a `GInteractions` object",
                    x=glue("`class(loops)` is {class(loops)}"),
                    i=glue("See `mariner::as_ginteractions()` to convert",
                           " a data-frame-like object to",
                           " a GInteractions object")))
        } else {
            rlang::abort(c("`loops` must be a `GInteractions` object",
                    x=glue("`class(loops)` is {class(loops)}")))
        }
    }

    ## Set scores to loops$score or abort if missing or wrong type/length
    if(missing(scores)){
        scores <- loops$score
        if(is.null(scores)){
            rlang::abort(c(glue("argument `scores` is missing, and `loops` ",
                         "does not contain a column named `score`")))
        }
    } else if(!is(scores, "numeric")){
        rlang::abort(c("`scores` must be a numeric vector",
                x = glue("`class(scores)` is {class(scores)}")))
    } else if(length(scores) != length(loops)){
        rlang::abort(c("`scores` must be the same length as `loops`",
                x = glue("`length(scores)` is {length(scores)} and ",
                         "`length(loops)` is {length(loops)}")))
    }

    ## Check that clusterType is valid
    possibleClusterTypes <- c("fast_greedy","walktrap", "leiden", "infomap",
                              "label_prop", "edge_betweenness",
                              "connected_components")
    if(!clusterType %in% possibleClusterTypes){
        rlang::abort(c(glue("`clusterType`=\"{clusterType}\" is not a valid",
                     " clustering algorithm option"),
                i=glue("`clusterType` must be one of \"",
                       glue_collapse({possibleClusterTypes}, sep ="\", \""),
                       "\"")))
    }

    ## Check that leiden resolution is a numeric
    if(!is(leidenResolution, "numeric")){
        rlang::abort(c("`leidenResolution` must be numeric",
                x=glue("`class(leidenResolution)` is ",
                       "{class(leidenResolution)}")))
    } else if(length(leidenResolution) > 1){
        rlang::warn(glue("More than one `leidenResolution` provided, only",
                         " the first will be used"))
        leidenResolution <- leidenResolution[1]
    }

    ## Give a warning if leidenResolution is provided but leiden is not
    ## used as clustering algorithm
    if(clusterType != "leiden"){
        if("leidenResolution" %in% names(call_args)){
            rlang::warn(glue("`leidenResolution` is provided but will not be",
                             " used because `clusterType` is not \"leiden\""))
        }
    }

    #### Filter Scores -------------------------------------------------------------
    ## calculate pruneUnder value if both original and added loops are present
    ## in "source" column
    if(missing(pruneUnder)){
        if(!is.null(loops$source)){
            tempLoops <- data.frame(source = loops$source,
                                    score = scores)
            subset <- split(tempLoops, tempLoops$source)
            originalLoops <- subset$original
            newLoops <- subset$added

            if(length(originalLoops) == 0 & length(newLoops == 0)){
                ## error if no value provided and source labels are unexpected
                rlang::abort(c(glue("`pruneUnder` could not be calculated, ",
                               "expected \"original\" and \"added\" ",
                               "in `source` column"),
                               i=glue("Use a `pruneUnder` value of 0 ",
                                "to keep all loops.")))
            }

            ## set pruneUnder to the median score of original loops
            pruneUnder = median(originalLoops$score)

        } else {
            ## error if no value provided and no source column
            rlang::abort(c("`pruneUnder` must be provided since there is no
                `source` column in `loops`",
                    i=glue("Use a `pruneUnder` value of 0 to keep all loops.")))
        }
    }

#### Build network -------------------------------------------------------------
    ## Remove scores from any added loops with a score less than pruneUnder
    if(!is.null(loops$source)){
        scores[scores < pruneUnder &
                   loops$source == "added"] <- NA
    } else {
        scores[scores < pruneUnder] <- NA
    }

    message(paste0("Pruning added loops with a score less than ",
                   pruneUnder))

    ## Set nodes to anchors and edges to scores
    anc <- InteractionSet::anchors(loops, type="both", id=TRUE)

    relations = data.frame(from = anc$first,
                           to = anc$second,
                           weights = scores,
                           chr = mariner::seqnames1(loops))

    relations[is.na(relations)] <- 0

    ## create community column, set to NA
    InteractionSet::regions(loops)$anchorCommunity <- NA

    ## keep track of total communities since community numbers will reset for
    ## each chromosome
    totalCommunities <- 0

    ## Create one graph per chrom
    chroms <- loops |>
        InteractionSet::regions() |>
        GenomicRanges::seqnames() |>
        unique()

    for(chrom in chroms){
        # filter to loops in chromosome with scores above pruneUnder
        relations_chrom <- dplyr::filter(relations,
                                         chr == chrom)
        g <- igraph::graph_from_data_frame(relations_chrom,
                                           directed=FALSE)

    #### Assign anchors to communities------------------------------------------
        # transform weights to log2
        relations_chrom$weights = log2(relations_chrom$weights)

        # remove edges with scores 0 or less
        g_filter <- igraph::delete_edges(g, which(relations_chrom$weights <= 0))
        relations_filter <- dplyr::filter(relations_chrom, weights > 0)

        ## generate communities using chosen algorithm
        if(clusterType == "edge_betweenness"){
            maxScore <- max(relations_filter$weights)
            flippedWeights <- 1 - ((relations_filter$weights)/maxScore)
            flippedWeights[which(flippedWeights==0)] <- 1e-10
        }

        between_g <-
            switch(clusterType,
                   fast_greedy =
                       igraph::cluster_fast_greedy(g_filter,
                                                   weights =
                                                       relations_filter$weights),
                   walktrap =
                       igraph::cluster_walktrap(g_filter,
                                                weights =
                                                    relations_filter$weights),
                   leiden =
                       igraph::cluster_leiden(g_filter,
                                              weights = relations_filter$weights,
                                              objective_function = "CPM",
                                              resolution_parameter =
                                                  leidenResolution,
                                              n_iterations = 5),
                   infomap =
                       igraph::cluster_infomap(g_filter,
                                               e.weights =
                                                   relations_filter$weights),
                   label_prop =
                       igraph::cluster_label_prop(g_filter,
                                                  weights =
                                                      relations_filter$weights),
                   edge_betweenness =
                       igraph::cluster_edge_betweenness(g_filter,
                                                        weights = flippedWeights,
                                                        directed = F),
                   connected_components =
                       igraph::components(g_filter)
            )

        ## put community numbers in genomic order
        anchorCommunities <-
            match(as.numeric(igraph::membership(between_g)),
              unique(as.numeric(igraph::membership(between_g))))

        ## assign loop anchors to communities
        # add totalCommunities to assigned communities to ensure unique
        # community numbers
        InteractionSet::regions(loops)$anchorCommunity[
            as.numeric(names(igraph::membership(between_g)))] <-
            anchorCommunities + totalCommunities

        ## add number of new unique communities to total communities
        totalCommunities <- totalCommunities +
            length(unique(as.numeric(igraph::membership(between_g))))
    }

#### adjust border nodes--------------------------------------------------------
    ## check if nodes on borders of communities fit well in neighboring
    ## communities based on scores
    anchors <- regions(loops) |> as.data.frame()

    bordersLeft <- anchors |>
        dplyr::group_by(anchorCommunity) |>
        dplyr::slice_head(n = 1)

    bordersRight <- anchors |>
        dplyr::group_by(anchorCommunity) |>
        dplyr::slice_tail(n = 1)

    ## TODO replace for loop with map
    communityNums <- unique(loops$loopCommunity) |>
        na.omit() |>
        as.numeric()

    movement <- GenomicRanges::GRanges()
    for(comm in communityNums){
        nodeLeft <- bordersLeft |>
            dplyr::filter(anchorCommunity == comm) |>
            plyranges::as_granges()

        nodeRight <- bordersRight |>
            dplyr::filter(anchorCommunity == comm) |>
            plyranges::as_granges()

        ln <- .compareNodeToComm(loops, nodeLeft, T)
        rn <- .compareNodeToComm(loops, nodeRight, F)

        movement <- c(movement, ln, rn)
    }

    ## filter to only anchors that moved
    if(length(movement) > 0){
        movement <- movement[which(sapply(movement$anchorCommunity,
                                          function(x) length(x) > 1))]
    }

    ## reassign anchor communities
    overlaps <- InteractionSet::findOverlaps(movement,
                                             InteractionSet::regions(loops))
    InteractionSet::regions(loops)[S4Vectors::subjectHits(overlaps)] <-
        movement[S4Vectors::queryHits(overlaps)]

    ## set loop community to any communities shared by both anchors
    anchor1com <- InteractionSet::regions(loops)[
        InteractionSet::anchorIds(loops)$first]$anchorCommunity
    anchor2com <- InteractionSet::regions(loops)[
        InteractionSet::anchorIds(loops)$second]$anchorCommunity

    loops$loopCommunity <-
        mapply(function(x, y) intersect(x,y), anchor1com, anchor2com)

    ## remove added loops under the pruning value
    loops <- loops[which(loops$source == "original" |
                             as.logical(scores >= pruneUnder))]

    S4Vectors::metadata(loops)$pruningValue <- pruneUnder

    return(loops)
}
