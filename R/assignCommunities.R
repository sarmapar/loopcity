#' Assign communities to interactions
#'
#' @param loops GInteractions object of interactions to include in the network
#' @param score a numeric vector of scores to use as edge weights,
#'  if `loops` contains a column named 'score', this column will be used
#' @param clusterType character vector of length 1 containing a keyword
#'  for the clustering algorithm to use. Must be one of "fast_greedy",
#'  "walktrap", "leiden", "infomap", "label_prop", or "edge_betweenness"
#' @return GInteractions object with community metadata
#' @export
#'
#' @examples
#' # Add example
assignCommunities <- function(loops,
                              score = loops$score,
                              clusterType = "leiden",
                              leiden_resolution = 0.01) {

  ## make network
  anc <- InteractionSet::anchors(loops, type="both", id=TRUE)

  relations = data.frame(from = anc$first,
                         to = anc$second,
                         weights = loops$score)
  relations[is.na(relations)] <- 0

  g <- igraph::graph_from_data_frame(relations,
                                directed=FALSE,
                                vertices=unique(c(anc[[2]],anc[[1]])))

  ## assign anchors to communities
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
           fast_greedy = igraph::cluster_fast_greedy(g_filter,
                                             weights = relations_filter$weights),
           walktrap = igraph::cluster_walktrap(g_filter,
                                       weights = relations_filter$weights),
           leiden = igraph::cluster_leiden(g_filter,
                                   weights = relations_filter$weights,
                                   objective_function = "CPM",
                                   resolution_parameter = leiden_resolution,
                                   n_iterations = 5),
           infomap = igraph::cluster_infomap(g_filter,
                                     e.weights = relations_filter$weights),
           label_prop = igraph::cluster_label_prop(g_filter,
                                          weights = relations_filter$weights),
           edge_betweenness = igraph::cluster_edge_betweenness(g_filter,
                                                   weights = flippedWeights,
                                                   directed = F)
  )

  InteractionSet::regions(loops)$anchorCommunity[as.numeric(names(igraph::membership(between_g)))] <-
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
