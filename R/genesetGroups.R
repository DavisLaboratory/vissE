#' Identify gene-set clusters from a gene-set overlap network
#'
#' This function identifies gene-set clusters from a gene-set overlap network
#' produced using vissE. Various graph clustering algorithms from the `igraph`
#' package can be used for clustering. Gene-set clusters identified are then
#' sorted based on their size and a given statistic of interest (absolute of the
#' statistic is maximised per cluster).
#'
#' @param ig an igraph object, containing a network of gene set overlaps computed
#'   using [computeMsigNetwork()].
#' @param alg a function, from the `igraph` package that should be used to
#'   perform graph-clustering (default is `igraph::cluster_walktrap`). The
#'   function should produce a `communities` object.
#' @param genesetStat a named numeric, containing statistics for each gene-set
#'   that are to be used in cluster prioritisation. If NULL, clusters are
#'   prioritised based on their size (number of gene-sets in them).
#' @param minSize a numeric, stating the minimum size a cluster can be (default
#'   is 2).
#' @param algparams a list, specifying additional parameters that are to be
#'   passed to the graph clustering algorithm.
#'
#' @return a list, containing gene-sets that belong to each cluster. Items in
#'   the list are organised based on prioritisation.
#' @export
#'
#' @details Gene-sets clusters are identified using graph clustering and are
#'   prioritised based on a combination of cluster size and optionally, a
#'   statistic of interest (e.g., enrichment scores). A product-of-ranks
#'   approach is used to prioritise clusters when gene-set statistics are
#'   available. In this approach, clusters are ranked based on their cluster
#'   size (largest to smallest) and on the median absolute statistic of
#'   gene-sets within it (largest to smallest). The product of these ranks is
#'   computed and clusters are ranked based on these product-of-rank statistic
#'   (smallest to largest).
#'
#'   When prioritising using cluster size and gene-set statistics, if statistics
#'   for some gene-sets in the network are missing, only the size is used in
#'   cluster prioritisation.
#'
#' @examples
#' data(hgsc)
#' ovlap <- computeMsigOverlap(hgsc, thresh = 0.25)
#' ig <- computeMsigNetwork(ovlap, hgsc)
#' findMsigClusters(ig)
findMsigClusters <- function(ig, genesetStat = NULL, minSize = 2, alg = igraph::cluster_walktrap, algparams = list()) {
  #param checks
  checkGraph(ig)
  checkNumericRange(minSize, 'minSize', pmin = 0)
  if (!is.null(genesetStat)) checkGenesetStat(genesetStat)
  
  #identify clusters
  emsg = 'graph clustering algorithm provided is invalid'
  algparams = c(list('graph' = ig), algparams)
  grps = tryCatch(do.call(alg, algparams), error = function(e) stop(emsg))
  if (!is(grps, 'communities')) {
    stop(emsg)
  }
  if (length(grps) == 0) {
    return(list())
  }
  
  #extract clustering results
  grps = igraph::groups(grps)
  
  #filter small groups
  grp.size = sapply(grps, length)
  grps = grps[grp.size >= minSize]
  grp.size = grp.size[grp.size >= minSize]
  
  #rank by cluster size
  rnk = rank(grp.size)
  if (!is.null(genesetStat)) {
    if (all(igraph::V(ig)$name %in% names(genesetStat))){
      #rank by cluster statistics
      rnk.stat = rank(sapply(grps, function(x) stats::median(abs(genesetStat[x]))))
      #combine and compute rank of product-of-ranks
      rnk = rank(rnk.stat * rnk)
    } else {
      warning("not using 'genesetStat' because statistics for some genesets are missing")
    }
  }
  grps = grps[order(rnk, decreasing = TRUE)]
  names(grps) = 1:length(grps)
  
  return(grps)
}