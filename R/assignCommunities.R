#' Assign communities to interactions
#'
#' @param loops GInteractions object of interactions to include in the network
#' @param scores a numeric vector of scores to use as edge weights,
#'  if `loops` contains a column named 'score', this column will be used
#' @param clusterType character vector of length 1 containing a keyword
#'  for the clustering algorithm to use. Must be one of "fast_greedy",
#'  "walktrap", "leiden", "infomap", "label_prop", or "edge_betweenness"
#' @param pruneUnder
#' @return GInteractions object with community metadata
#'
#' @importFrom rlang abort
#' @importFrom glue glue glue_collapse
#' @export
#'
#' @examples
#' # Add example
assignCommunities <- function(loops,
                              scores,
                              clusterType = "leiden",
                              leiden_resolution = 0.01,
                              pruneUnder) {

## Parameter checking ----------------------------------------------------------
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
                              "label_prop", "edge_betweenness")
    if(!clusterType %in% possibleClusterTypes){
        rlang::abort(c(glue("`clusterType`=\"{clusterType}\" is not a valid",
                     " clustering algorithm option"),
                i=glue("`clusterType` must be one of \"",
                       glue_collapse({possibleClusterTypes}, sep ="\", \""),
                       "\"")))
    }

    ## Check that leiden resolution is a numeric
    if(!is(leiden_resolution, "numeric")){
        rlang::abort(c("`leiden_resolution` must be numeric",
                x=glue("`class(leiden_resolution)` is ",
                       "{class(leiden_resolution)}")))
    } else if(length(leiden_resolution) > 1){
        rlang::warn(glue("More than one `leiden_resolution` provided, only",
                         " the first will be used"))
        leiden_resolution <- leiden_resolution[1]
    }

    ## Give a warning if leiden_resolution is provided but leiden is not
    ## used as clustering algorithm
    if(clusterType != "leiden"){
        if(leiden_resolution != 0.01){
            rlang::warn(glue("`leiden_resolution` is provided but will not be",
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

            lowerLimit <- quantile(newLoops$score, 0.4)
            upperLimit <- quantile(originalLoops$score, 0.6)
            ogLoopDensity <- density(as.numeric(originalLoops$score),
                                     from = lowerLimit,
                                     to = upperLimit,
                                     n = 2^10)
            newLoopDensity <- density(as.numeric(newLoops$score),
                                      from = lowerLimit,
                                      to = upperLimit,
                                      n = 2^10)

            densityDiff <- ogLoopDensity$y - newLoopDensity$y
            intersectionPoint <- newLoopDensity$x[which(
                diff(densityDiff > 0) != 0) + 1]

            if(length(intersectionPoint) == 1){
                pruneUnder <- intersectionPoint
            } else{
                rlang::warn(glue("Optimal pruning threshold could not be found",
                                 "setting `pruneUnder` to 1."))
                pruneUnder <- 1
            }
            message(paste0("Pruning all loops with a score less than ", pruneUnder))

        } else {
            ## error if no value provided and no added loops
            rlang::abort(c("`pruneUnder` must be provided since there is no
                `source` column in `loops`",
                    i=glue("Use a `pruneUnder` value of 0 to keep all loops.")))
        }
    }

#### Build network -------------------------------------------------------------
    anc <- InteractionSet::anchors(loops, type="both", id=TRUE)

    relations = data.frame(from = anc$first,
                           to = anc$second,
                           weights = scores)
    relations[is.na(relations)] <- 0

    g <- igraph::graph_from_data_frame(relations,
                                       directed=FALSE,
                                       vertices=unique(c(anc[[2]],anc[[1]])))

#### assign anchors to communities----------------------------------------------
    # remove edges with scores less than pruneUnder
    g_filter <- igraph::delete_edges(g, which(relations$weights <= pruneUnder))
    relations_filter <- dplyr::filter(relations, weights > pruneUnder)

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
                                              leiden_resolution,
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
                                                    directed = F)
        )

#### adjust border nodes--------------------------------------------------------
    InteractionSet::regions(loops)$anchorCommunity[
        as.numeric(names(igraph::membership(between_g)))] <-
        as.numeric(igraph::membership(between_g))

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

    ## re-prune loops under the pruning value
    loops$loopCommunity <- ifelse(loops$score >= pruneUnder,
                                  loops$loopCommunity,
                                  list(numeric(0)))

    metadata(loops)$pruningValue <- pruneUnder

    return(loops)
}
