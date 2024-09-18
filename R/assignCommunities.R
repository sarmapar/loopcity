#' Assign communities to interactions
#'
#' @param loops GInteractions object of interactions to include in the network
#' @param scores a numeric vector of scores to use as edge weights,
#'  if `loops` contains a column named 'score', this column will be used
#' @param clusterType character vector of length 1 containing a keyword
#'  for the clustering algorithm to use. Must be one of "fast_greedy",
#'  "walktrap", "leiden", "infomap", "label_prop", or "edge_betweenness"
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
                              leiden_resolution = 0.01) {

## Parameter checking ----------------------------------------------------------
    ## Check if loops is a GInteractions object, with helpful message for
    ## converting data frame to GInteractions
    if(!is(loops,"GInteractions")){
        if(is(loops,"data.frame") | is(loops,"data.table") |
                                       is(loops,"DataFrame")){
            abort(c("`loops` must be a `GInteractions` object",
                    x=glue("`class(loops)` is {class(loops)}"),
                    i=glue("See `mariner::as_ginteractions()` to convert",
                           " a data-frame-like object to",
                           " a GInteractions object")))
        } else {
            abort(c("`loops` must be a `GInteractions` object",
                    x=glue("`class(loops)` is {class(loops)}")))
        }
    }

    ## Set scores to loops$score or abort if missing or wrong type/length
    if(missing(scores)){
        scores <- loops$score
        if(is.null(scores)){
            abort(c(glue("argument `scores` is missing, and `loops` ",
                         "does not contain a column named `score`")))
        }
    } else if(!is(scores, "numeric")){
        abort(c("`scores` must be a numeric vector",
                x = glue("`class(scores)` is {class(scores)}")))
    } else if(length(scores) != length(loops)){
        abort(c("`scores` must be the same length as `loops`",
                x = glue("`length(scores)` is {length(scores)} and ",
                         "`length(loops)` is {length(loops)}")))
    }

    ## Check that clusterType is valid
    possibleClusterTypes <- c("fast_greedy","walktrap", "leiden", "infomap",
                              "label_prop", "edge_betweenness")
    if(!clusterType %in% possibleClusterTypes){
        abort(c(glue("`clusterType`=\"{clusterType}\" is not a valid",
                     " clustering algorithm option"),
                i=glue("`clusterType` must be one of \"",
                       glue_collapse({possibleClusterTypes}, sep ="\", \""),
                       "\"")))
    }

    ## Check that leiden resolution is a numeric
    if(!is(leiden_resolution, "numeric")){
        abort(c("`leiden_resolution` must be numeric",
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
    # remove edges with 0
    g_filter <- igraph::delete_edges(g, which(relations$weights <= 0))
    relations_filter <- dplyr::filter(relations, weights > 0)

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

    InteractionSet::regions(loops)$anchorCommunity[
        as.numeric(names(igraph::membership(between_g)))] <-
        as.numeric(igraph::membership(between_g))

    ## if both anchors of a loop are in the same community,
    ## set the loop community
    ## otherwise set loopCommunity to NA
    anchor1com <- InteractionSet::regions(loops)[
        InteractionSet::anchorIds(loops)$first]$anchorCommunity
    anchor2com <- InteractionSet::regions(loops)[
        InteractionSet::anchorIds(loops)$second]$anchorCommunity

    loops$loopCommunity <-
        ifelse(anchor1com == anchor2com,
               anchor1com,
               NA)

    return(loops)
}
