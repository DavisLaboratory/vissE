#' @importFrom igraph V
NULL

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
characteriseGeneset <- function(gs, thresh = 0.15, measure = c('jaccard', 'ovlapcoef')) {
  measure = match.arg(measure)
  gsc = msigdb::msigdb.hs.SYM()
  gsc = msigdb::appendKEGG(gsc)
  
  #filter out large and small gene sets
  len = sapply(lapply(gsc, GSEABase::geneIds), length)
  gsc = GSEABase::GeneSetCollection(gsc[len > 10 & len < 500])

  #compute overlaps
  ovmat = computeMsigOverlap(gsc, GSEABase::GeneSetCollection(gs), 0.3, 'ovlapcoef')
  gsc = GSEABase::GeneSetCollection(c(gs, gsc))
  
  #identify neighbours
  ovmat = ovmat[ovmat$weight > thresh, ]
  
  #induce graph
  nb = GSEABase::GeneSetCollection(gsc[nb])
  ovmat = computeMsigOverlap(nb, thresh, measure)
  nbnet = computeMsigNetwork(ovmat, gsc)
  
  return(nbnet)
}
