#' @importFrom igraph E V E<- V<-
NULL

#' Compute gene set overlap
#'
#' Compute overlap between gene sets from a GeneSetCollection using the Jaccard
#' index or the overlap coefficient. These values can then be used to compute a
#' network of gene set overlaps.
#'
#' @param msigGsc1 a GeneSetCollection object.
#' @param msigGsc2 a GeneSetCollection object or NULL if pairwise overlaps are
#'   to be computed.
#' @param thresh a numeric, specifying the threshold to discard pairs of gene
#'   sets.
#' @param measure a character, specifying the similarity measure to use:
#'   `jaccard` for the Jaccard Index and `ovlapcoef` for the Overlap
#'   Coefficient.
#'
#' @return a data.frame, containing the overlap structure of gene sets
#'   represented as a network in the simple interaction format (SIF).
#' @export
#' @examples
#' data(hgsc)
#' ovlap <- computeMsigOverlap(hgsc)
#'
computeMsigOverlap <- function(msigGsc1, msigGsc2 = NULL, thresh = 0.25, measure = c('jaccard', 'ovlapcoef')) {
  #check collection size
  stopifnot(length(msigGsc1) > 0)
  stopifnot(all(sapply(lapply(msigGsc1, GSEABase::geneIds), length) > 0))
  if (!is.null(msigGsc2)) {
    stopifnot(length(msigGsc2) > 0)
    stopifnot(all(sapply(lapply(msigGsc2, GSEABase::geneIds), length) > 0))
  }
  #check threshold
  stopifnot(thresh >= 0 & thresh <= 1)
  measure = match.arg(measure)
  
  #combine genesets
  gsc = msigGsc1
  is1 = rep(TRUE, length(msigGsc1))
  if (!is.null(msigGsc2)) {
    #get sets unique to msigGsc2
    gsc2 = msigGsc2[setdiff(names(msigGsc2), names(gsc))]
    gsc = GSEABase::GeneSetCollection(c(gsc, gsc2))
    is2 = names(gsc) %in% names(msigGsc2)
  }
  
  #compute incidence matrix and split
  imat = GSEABase::incidence(gsc)
  if (is.null(msigGsc2)) {
    ovlap = tcrossprod(imat)
    len1 = len2 = rowSums(imat)
  } else {
    ovlap = tcrossprod(imat[is1, , drop = FALSE], imat[is2, , drop = FALSE])
    len = rowSums(imat)
    len1 = len[is1]
    len2 = len[is2]
  }
  
  #overlap coef
  if (measure %in% 'jaccard') {
    mat = ovlap / (outer(len1, len2, '+') - ovlap)
  } else {
    mat = ovlap / outer(len1, len2, pmin)
  }

  #convert to data.frame
  mat = reshape2::melt(mat, varnames = c('gs1', 'gs2'), value.name = 'weight')
  mat$gs1 = as.character(mat$gs1)
  mat$gs2 = as.character(mat$gs2)
  mat = mat[mat$weight >= thresh, , drop = FALSE]

  if (is.null(msigGsc2)) {
    #remove symmetric values
    mat = mat[mat$gs1 < mat$gs2, , drop = FALSE]
  }

  rownames(mat) = NULL

  return(mat)
}

#' Compute a network using computed gene set overlap
#'
#' Computes an igraph object using information on gene sets and gene sets
#' computed using the [computeMsigOverlap()] function.
#'
#' @param genesetOverlap a data.frame, containing results of an overlap analysis
#'   computed using the [computeMsigOverlap()] function.
#' @param msigGsc a GeneSetCollection object, containing gene sets used to
#'   compute overlap.
#'
#' @return an igraph object
#' @export
#'
#' @examples
#' data(hgsc)
#' ovlap <- computeMsigOverlap(hgsc)
#' ig <- computeMsigNetwork(ovlap, hgsc)
#'
computeMsigNetwork <- function(genesetOverlap, msigGsc) {
  stopifnot(nrow(genesetOverlap) > 0)

  #select genesets in the network
  setnames = sapply(msigGsc, GSEABase::setName)
  setnames = intersect(setnames, c(genesetOverlap$gs1, genesetOverlap$gs2))
  stopifnot(all(genesetOverlap$gs1 %in% setnames))
  stopifnot(all(genesetOverlap$gs2 %in% setnames))
  msigGsc = msigGsc[setnames]

  #check Broad collection
  isBroad = sapply(msigGsc, function(x)
    GSEABase::collectionType(GSEABase::collectionType(x)) %in% 'Broad')

  #compute nodedf
  ovNodes = data.frame(
    'Name' = setnames,
    'Size' = sapply(lapply(msigGsc, GSEABase::geneIds), length),
    'Category' = 'custom',
    'SubCategory' = NA
  )
  ovNodes$Category[isBroad] = unlist(sapply(sapply(msigGsc[isBroad], GSEABase::collectionType),
                                     GSEABase::bcCategory))
  ovNodes$SubCategory[isBroad] = unlist(sapply(sapply(msigGsc[isBroad], GSEABase::collectionType),
                                        GSEABase::bcSubCategory))

  #create igraph
  msig_ig = igraph::graph_from_data_frame(genesetOverlap, directed = FALSE, vertices = ovNodes)

  return(msig_ig)
}

getMsigNeighbour <- function(srcsig, ig, thresh = 0.15) {
  stopifnot(thresh >= 0 & thresh <= 1)
  stopifnot(srcsig %in% V(ig)$name)

  #select surrounding network
  sub_ig = igraph::induced_subgraph(ig, vids = igraph::neighborhood(ig, nodes = srcsig)[[1]])
  if (all(E(sub_ig)$weight <= thresh))
    return(c())

  #extract confident edges
  sub_ig = igraph::subgraph.edges(sub_ig, E(sub_ig)[E(sub_ig)$weight > thresh])
  if (!srcsig %in% V(sub_ig)$name)
    return(c())

  return(igraph::neighbors(sub_ig, srcsig)$name)
}
