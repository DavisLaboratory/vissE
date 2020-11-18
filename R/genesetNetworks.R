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
#'   represented as a network in the simple interaction format (SIF)
#' @export
#' @examples
#'
#' data(hgsc)
#' ovlap <- computeMsigOverlap(hgsc)
#'
computeMsigOverlap <- function(msigGsc1, msigGsc2 = NULL, thresh = 0.1, measure = c('jaccard', 'ovlapcoef')) {
  stopifnot(thresh >= 0 & thresh <= 1)
  measure = match.arg(measure)

  #return empty result if no gene sets are provided
  if (length(msigGsc1) == 0 | !is.null(msigGsc2) & length(msigGsc2) == 0)
    return(data.frame('gs1' = character(), 'gs2' = character(), 'coef' = numeric()))

  #filter out very small and very large genesets
  genes1 = lapply(msigGsc1, GSEABase::geneIds)
  names(genes1) = sapply(msigGsc1, GSEABase::setName)
  len1 = sapply(genes1, length)
  genes1 = genes1[len1 > 10 & len1 < 500]
  len1 = len1[len1 > 10 & len1 < 500]

  #compute overlap network
  if (!is.null(msigGsc2)) {
    #filter out very small and very large genesets
    genes2 = lapply(msigGsc2, GSEABase::geneIds)
    names(genes2) = sapply(msigGsc2, GSEABase::setName)
    len2 = sapply(genes2, length)
    genes2 = genes2[len2 > 10 & len2 < 500]
    len2 = len2[len2 > 10 & len2 < 500]
  } else {
    genes2 = NULL
    len2 = len1
  }
  ovlap = intersectSize(genes1, genes2)

  #overlap coef
  if (measure %in% 'jaccard') {
    mat = ovlap / (outer(len1, len2, '+') - ovlap)
  } else {
    mat = ovlap / outer(len1, len2, pmin)
  }

  #convert to data.frame
  mat = reshape2::melt(mat, varnames = c('gs1', 'gs2'), value.name = 'coef')
  mat$gs1 = as.character(mat$gs1)
  mat$gs2 = as.character(mat$gs2)
  mat = mat[mat$gs1 < mat$gs2 & mat$coef >= thresh, , drop = FALSE]
  rownames(mat) = NULL

  return(mat)
}

intersectSize <- function(x, y = NULL) {
  vals = unique(unlist(c(x, y)))

  #compute overlap network
  matx = plyr::laply(x, function(s) as.numeric(vals %in% s))
  rownames(matx) = names(x)
  colnames(matx) = vals

  if (is.null(y)) {
    ovlap = tcrossprod(matx)
  } else {
    maty = plyr::laply(y, function(s) as.numeric(vals %in% s))
    rownames(maty) = names(y)
    colnames(maty) = vals
    ovlap = tcrossprod(matx, maty)
  }

  return(ovlap)
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
#'
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

#' Get similar gene sets from overlap networks
#'
#' @param srcsig a character, naming the gene set to begin the search with.
#' @param ig an igraph object, containing a network of gene set overlaps
#'   computed using [computeMsigNetwork()].
#' @param thresh a numeric, specifying the threshold to discard pairs of gene
#'   sets.
#'
#' @return a character, containing the names of gene sets that overlap with the
#'   source signature.
#' @export
#'
#' @examples
#'
#' data("msigOverlapNetwork")
#' neighbours <- getMsigNeighbour('HALLMARK_HYPOXIA', msigOverlapNetwork, 0.1)
#'
getMsigNeighbour <- function(srcsig, ig, thresh = 0) {
  stopifnot(thresh >= 0 & thresh <= 1)
  stopifnot(srcsig %in% V(ig)$name)

  #select surrounding network
  sub_ig = igraph::induced_subgraph(ig, vids = igraph::neighborhood(ig, nodes = srcsig)[[1]])
  if (all(E(sub_ig)$coef <= thresh))
    return(c())

  #extract confident edges
  sub_ig = igraph::subgraph.edges(sub_ig, E(sub_ig)[E(sub_ig)$coef > thresh])
  if (!srcsig %in% V(sub_ig)$name)
    return(c())

  return(igraph::neighbors(sub_ig, srcsig)$name)
}
