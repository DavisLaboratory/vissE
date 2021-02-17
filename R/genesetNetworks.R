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
computeMsigOverlap <- function(msigGsc1, msigGsc2 = NULL, thresh = 0.15, measure = c('jaccard', 'ovlapcoef')) {
  stopifnot(thresh >= 0 & thresh <= 1)
  measure = match.arg(measure)
  
  #empty genesets
  gsc1_lengths = sapply(lapply(msigGsc1, GSEABase::geneIds), length)
  stopifnot(all(gsc1_lengths > 0))
  
  #get organism
  idType1 = msigdb::getMsigIdType(msigGsc1)
  organism1 = msigdb::getMsigOrganism(msigGsc1, idType1)

  if (!is.null(msigGsc2)) {
    #empty genesets
    gsc2_lengths = sapply(lapply(msigGsc2, GSEABase::geneIds), length)
    stopifnot(all(gsc2_lengths > 0))
    
    #get organism
    idType2 = msigdb::getMsigIdType(msigGsc2)
    organism2 = msigdb::getMsigOrganism(msigGsc1, idType2)
    
    #check for concordance
    stopifnot(idType1 == idType2)
    stopifnot(organism1 == organism2)
  }

  #return empty result if no gene sets are provided
  if (length(msigGsc1) == 0 | !is.null(msigGsc2) & length(msigGsc2) == 0)
    return(data.frame('gs1' = character(), 'gs2' = character(), 'weight' = numeric()))

  #filter out very small and very large genesets
  genes1 = lapply(msigGsc1, GSEABase::geneIds)
  names(genes1) = sapply(msigGsc1, GSEABase::setName)
  len1 = sapply(genes1, length)
  # genes1 = genes1[len1 > 10 & len1 < 500]
  # len1 = len1[len1 > 10 & len1 < 500]

  #compute overlap network
  if (!is.null(msigGsc2)) {
    #filter out very small and very large genesets
    genes2 = lapply(msigGsc2, GSEABase::geneIds)
    names(genes2) = sapply(msigGsc2, GSEABase::setName)
    len2 = sapply(genes2, length)
    # genes2 = genes2[len2 > 10 & len2 < 500]
    # len2 = len2[len2 > 10 & len2 < 500]
  } else {
    genes2 = NULL
    len2 = len1
  }
  ovlap = intersectSize(genes1, genes2, organism1, idType1)
  ovlap = as.matrix(ovlap)

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

intersectSize <- function(x, y = NULL, org, idType) {
  allg = unique(unlist(c(x, y)))

  #retrieve/create membership matrix for x
  matx = retrieveMat(x, allg, org, idType)

  if (is.null(y)) {
    #compute overlap
    ovlap = Matrix::tcrossprod(matx)
  } else {
    #retreive/create membership matrix for x
    maty = retrieveMat(y, allg, org, idType)
    
    #compute overlap
    ovlap = Matrix::tcrossprod(matx, maty)
  }

  return(ovlap)
}

retrieveMat <- function(gslist, allg, org, idType) {
  if (org %in% 'hs') {
    mem_mat = mem_mat_hs
  } else {
    mem_mat = mem_mat_mm
  }
  
  if (class(idType) %in% 'SymbolIdentifier') {
    if (org %in% 'hs') {
      colnames(mem_mat) = as.character(
        AnnotationDbi::mapIds(
          org.Hs.eg.db::org.Hs.eg.db,
          colnames(mem_mat),
          'SYMBOL',
          'ENTREZID'
        )
      )
    } else {
      colnames(mem_mat) = as.character(
        AnnotationDbi::mapIds(
          org.Mm.eg.db::org.Mm.eg.db,
          colnames(mem_mat),
          'SYMBOL',
          'ENTREZID'
        )
      )
    }
  }
  
  #retrieve precomputed results
  mem_mat = mem_mat[intersect(names(gslist), rownames(mem_mat)), intersect(allg, colnames(mem_mat))]
  
  #determine new genes and gene sets
  newgs = setdiff(names(gslist), rownames(mem_mat))
  allg = setdiff(allg, colnames(mem_mat))
  
  #compute overlap for new gene sets
  if (length(newgs) > 0) {
    #initialise
    gids = colnames(mem_mat)
    mat = Matrix::Matrix(
      0,
      nrow = length(newgs),
      ncol = ncol(mem_mat),
      dimnames = list(newgs, gids),
      sparse = TRUE
    )
    #compute membership
    for (i in 1:length(newgs)) {
      mat[i, ] = as.numeric(gids %in% gslist[[i]])
    }
    #merge
    mem_mat = rbind(mem_mat, mat)
  }
  
  #compute overlap for new genes
  if (length(allg) > 0) {
    #initialise
    gsnames = rownames(mem_mat)
    mat = Matrix::Matrix(
      0,
      nrow = length(gsnames),
      ncol = length(allg),
      dimnames = list(gsnames, allg),
      sparse = TRUE
    )
    #compute membership
    for (i in 1:length(gsnames)) {
      mat[i, ] = as.numeric(allg %in% gslist[[i]])
    }
    #merge
    mem_mat = cbind(mem_mat, mat)
  }
  
  return(mem_mat)
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
