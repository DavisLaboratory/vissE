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
#' @param gscolcs a character, listing the MSigDB collections to use as a
#'   background (defaults to h, c2, and c5). Collection types can be retrieved
#'   using [msigdb::listCollections()].
#'
#' @inheritParams computeMsigOverlap
#'
#' @return an igraph object, containing gene sets that are similar to the query
#'   set. The network contains relationships between results of the query too.
#' @export
#'
#' @examples
#' library(GSEABase)
#' data(hgsc)
#'
#' #create a geneset using one of the Hallmark gene sets
#' mySet <- GeneSet(
#'   geneIds(hgsc[[2]]),
#'   setName = 'MySet',
#'   geneIdType = SymbolIdentifier()
#' )
#'
#' \donttest{
#' #characterise the custom gene set
#' ig <- characteriseGeneset(mySet)
#' plotMsigNetwork(ig)
#' }
#' 
characteriseGeneset <- function(gs, thresh = 0.2, measure = c('ovlapcoef', 'jaccard'), gscolcs = c('h', 'c2', 'c5')) {
  measure = match.arg(measure)
  
  #retrieve appropriate GeneSetCollection
  gsc_gs = GSEABase::GeneSetCollection(gs)
  id = msigdb::getMsigIdType(gsc_gs)
  org = msigdb::getMsigOrganism(gsc_gs, id)
  if (org %in% 'hs') {
    if (is(id, 'SymbolIdentifier')) {
      gsc = msigdb::msigdb.v7.2.hs.SYM()
    } else {
      gsc = msigdb::msigdb.v7.2.hs.EZID()
    }
  } else {
    if (is(id, 'SymbolIdentifier')) {
      gsc = msigdb::msigdb.v7.2.mm.SYM()
    } else {
      gsc = msigdb::msigdb.v7.2.mm.EZID()
    }
  }
  gsc = msigdb::appendKEGG(gsc)
  
  #subset collections to use
  stopifnot(all(gscolcs %in% msigdb::listCollections(gsc)))
  gsc = msigdb::subsetCollection(gsc, gscolcs)
  
  #filter out large and small gene sets
  len = sapply(lapply(gsc, GSEABase::geneIds), length)
  gsc = GSEABase::GeneSetCollection(gsc[len > 10 & len < 500])

  #compute overlaps
  ovmat = computeMsigOverlap(gsc, GSEABase::GeneSetCollection(gs), thresh, measure)
  gsc = GSEABase::GeneSetCollection(c(gs, gsc))
  
  #identify neighbours
  ovmat = ovmat[ovmat$weight > thresh, ]
  
  #induce graph
  nb = GSEABase::GeneSetCollection(gsc[ovmat[, 1]])
  ovmat = computeMsigOverlap(nb, thresh = 0.15, measure = 'jaccard')
  nbnet = computeMsigNetwork(ovmat, gsc)
  
  return(nbnet)
}
