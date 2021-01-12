#' Functionally characterise a list of genes
#'
#' This function can be used to perform a network-based enrichment analysis of a
#' list of genes. The list of genes are characterised based on their similarity
#' with gene sets from the MSigDB. A network of similar gene sets is retrieved
#' using this function.
#'
#' @param gs a GeneSet object, representing the list of genes that need to be
#'   characterised.
#'
#' @inheritParams computeMsigOverlap
#'
#' @return an igraph object, containing gene sets that are similar to the query
#'   set. The network contains relationships between results of the query too.
#' @export
#'
#' @examples
#'
#' data(hgsc)
#'
#' #create a geneset using one of the Hallmark gene sets
#' mySet <- GSEABase::GeneSet(GSEABase::geneIds(hgsc[[2]]), setName = 'MySet')
#'
#' \dontrun{
#' #characterise the custom gene set
#' ig <- characteriseGeneset(mySet)
#' plotMsigNetwork(ig)
#' }
#'
characteriseGeneset <- function(gs, thresh = 0.1, measure = c('ovlapcoef', 'jaccard')) {
  measure = match.arg(measure)

  #compute overlaps
  ovmat = computeMsigOverlap(GSEABase::GeneSetCollection(gs), msigdb, thresh, measure)
  gsc = GSEABase::GeneSetCollection(c(gs, msigdb))

  #combine new graph with precomputed graph
  nodedf = igraph::as_data_frame(msigOverlapNetwork, 'vertices')[, -(5:7)]
  edgedf = igraph::as_data_frame(msigOverlapNetwork, 'edges')
  colnames(ovmat) = colnames(edgedf)
  edgedf = rbind(edgedf, ovmat)

  newnodes = setdiff(union(edgedf$from, edgedf$to), c(nodedf$name, GSEABase::setName(gs)))
  newnodes = msigdb[newnodes]
  newnodedf = data.frame(
    sapply(c(newnodes, gs), GSEABase::setName),
    sapply(lapply(c(newnodes, gs), GSEABase::geneIds), length),
    c(sapply(lapply(newnodes, GSEABase::collectionType), GSEABase::bcCategory), 'custom'),
    c(sapply(lapply(newnodes, GSEABase::collectionType), GSEABase::bcSubCategory), 'custom')
  )
  colnames(newnodedf) = colnames(nodedf)
  nodedf = rbind(nodedf, newnodedf)

  fullig = igraph::graph_from_data_frame(edgedf, FALSE, nodedf)

  #define group by identifying communities
  nb = getMsigNeighbour(GSEABase::setName(gs), fullig, thresh)
  nbnet = igraph::induced_subgraph(fullig, vids = nb)

  return(nbnet)
}
