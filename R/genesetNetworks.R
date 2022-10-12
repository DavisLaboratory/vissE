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
#' @param measure a character, specifying the similarity measure to use: `ari`
#'   for the Adjusted Rand Index, `jaccard` for the Jaccard Index and
#'   `ovlapcoef` for the Overlap Coefficient.
#'
#' @return a data.frame, containing the overlap structure of gene sets
#'   represented as a network in the simple interaction format (SIF).
#' @export
#' @examples
#' data(hgsc)
#' ovlap <- computeMsigOverlap(hgsc)
#' 
computeMsigOverlap <- function(msigGsc1, msigGsc2 = NULL, thresh = 0.25, measure = c('ari', 'jaccard', 'ovlapcoef')) {
  #param checks
  checkGenesetCollection(msigGsc1)
  if (!is.null(msigGsc2)) checkGenesetCollection(msigGsc2)
  
  #check threshold
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
  total = ncol(imat)
  if (measure %in% 'ari') {
    mat = overlap.ari(len1, len2, ovlap, total)
  } else if (measure %in% 'jaccard') {
    mat = overlap.jaccard(len1, len2, ovlap, total)
  } else {
    mat = overlap.ovlapcoef(len1, len2, ovlap, total)
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

overlap.ari <-  function(len1, len2, ovlap, total) {
  #----design----
  # n = length(v1) # total genes
  # x = sum(v1 & v2) # intersection
  # a = sum(v1) # length of geneset 1
  # b = sum(v2) # length of geneset 2
  # 
  # sum_a = c(a, n - a)
  # sum_a = sum(sum_a * (sum_a - 1))
  # 
  # sum_b = c(b, n - b)
  # sum_b = sum(sum_b * (sum_b - 1))
  # 
  # prod_ab = sum_a * sum_b / (n * (n - 1))
  # 
  # sum_nij = c(x, n - (a + b - x), a - x, b - x)
  # sum_nij = sum(sum_nij * (sum_nij - 1))
  # 
  # (sum_nij - prod_ab) / ((sum_a + sum_b) / 2 - prod_ab)
  
  #define components of the equation
  n = total # total genes
  x = ovlap # intersection of gene-sets
  a = len1 # lengths of the first gene-set
  b = len2 # lengths of the second gene-set
  
  #compute sums of pairs within and outside each gene-set (first)
  sum_a = cbind(a, n - a) * cbind(a - 1, n - a - 1)
  sum_a = sum_a[, 1] + sum_a[, 2]
  sum_b = cbind(b, n - b) * cbind(b - 1, n - b - 1)
  sum_b = sum_b[, 1] + sum_b[, 2]
  sum_ab = outer(sum_a, sum_b, '+')
  prod_ab = outer(sum_a, sum_b, '*') / (n * (n - 1))
  
  #compute sums of shared pairs within or outside two gene-sets
  sum_nij = x * (x - 1)
  tmp = n - outer(a, b, '+') + x
  sum_nij = sum_nij + (tmp * (tmp - 1))
  tmp = a - x
  sum_nij = sum_nij + (tmp * (tmp - 1))
  tmp = t(b - t(x))
  sum_nij = sum_nij + (tmp * (tmp - 1))
  
  #ARI
  ari = (sum_nij - prod_ab) / ((sum_ab) / 2 - prod_ab)
  
  return(ari)
}

overlap.jaccard <- function(len1, len2, ovlap, total) {
  mat = ovlap / (outer(len1, len2, '+') - ovlap)
  return(mat)
}

overlap.ovlapcoef <- function(len1, len2, ovlap, total) {
  mat = ovlap / outer(len1, len2, pmin)
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
  #check params
  checkGenesetCollection(msigGsc)
  
  if (!nrow(genesetOverlap) > 0)
    stop("'genesetOverlap' should not be empty")

  #check for unknown gene-set names
  ovlapnames = c(genesetOverlap$gs1, genesetOverlap$gs2)
  unknownSets = setdiff(ovlapnames, sapply(msigGsc, GSEABase::setName))
  if (length(unknownSets > 0)) {
    unknownSets = paste(unknownSets, collapse = ', ')
    stop("the following gene-set names in 'genesetOverlap' are missing 'msigGsc': %s", unknownSets)
  }
  
  #select genesets in the network
  setnames = setdiff(ovlapnames, unknownSets)
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
  if(!srcsig %in% V(ig)$name)
    stop(sprintf("'%s' missing in the graph", srcsig))

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
