#' Compute gene set overlap
#'
#' Compute overlap between gene sets from a GeneSetCollection using the Jaccard
#' index or the overlap coefficient. These values can then be used to compute a
#' network of gene set overlaps.
#'
#' @param msigGsc a GeneSetCollection object.
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
computeMsigOverlap <- function(msigGsc, thresh = 0.1, measure = c('jaccard', 'ovlapcoef')) {
  stopifnot(thresh >= 0 & thresh <= 1)
  measure = match.arg(measure)

  #return empty result if no gene sets are provided
  if (length(msigGsc) == 0)
    return(data.frame('gs1' = character(), 'gs2' = character(), 'coef' = numeric()))

  #filter out very small and very large genesets
  genes = lapply(msigGsc, GSEABase::geneIds)
  names(genes) = sapply(msigGsc, GSEABase::setName)
  alllen = sapply(genes, length)
  genes = genes[alllen > 10 & alllen < 500]
  alllen = alllen[alllen > 10 & alllen < 500]

  #compute overlap network
  allg = unique(unlist(genes))
  gmat = plyr::laply(genes, function(x) as.numeric(allg %in% x))
  rownames(gmat) = names(genes)
  colnames(gmat) = allg
  ovlap = tcrossprod(gmat)

  #overlap coef
  if (measure %in% 'jaccard') {
    mat = ovlap / (outer(alllen, alllen, '+') - ovlap)
  } else {
    mat = ovlap / outer(alllen, alllen, pmin)
  }

  #convert to data.frame
  mat[!lower.tri(mat)] = NA
  mat = reshape2::melt(mat, varnames = c('gs1', 'gs2'), value.name = 'coef')
  mat = mat[!is.na(mat$coef), , drop = FALSE]
  mat = mat[mat$coef >= thresh, , drop = FALSE]
  rownames(mat) = NULL
  mat$gs1 = as.character(mat$gs1)
  mat$gs2 = as.character(mat$gs2)

  return(mat)
}

intersectSize <- function(x, y = NULL) {
  #compute overlap network
  vals = unique(unlist(c(x, y)))
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
  ovNodes$Category[isBroad] = sapply(sapply(msigGsc[isBroad], GSEABase::collectionType),
                                     GSEABase::bcCategory)
  ovNodes$SubCategory[isBroad] = sapply(sapply(msigGsc[isBroad], GSEABase::collectionType),
                                        GSEABase::bcSubCategory)

  #create igraph
  msig_ig = igraph::graph_from_data_frame(genesetOverlap, directed = FALSE, vertices = ovNodes)

  return(msig_ig)
}

