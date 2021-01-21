#' The Hallmark collection from the MSigDB
#'
#' The molecular signatures database (MSigDB) is a collection of over 25000 gene
#' expression signatures. Signatures in v7.2 are divided into 9 categories. The
#' Hallmarks collection contains gene expression signatures representing
#' molecular processes that are hallmarks in cancer development and progression.
#'
#' @format A GeneSetCollection object with 50 GeneSet objects representing the
#'   50 Hallmark gene expression signatures.
#' @docType data
#' @references Subramanian, A., Tamayo, P., Mootha, V. K., Mukherjee, S., Ebert,
#'   B. L., Gillette, M. A., ... & Mesirov, J. P. (2005). Gene set enrichment
#'   analysis: a knowledge-based approach for interpreting genome-wide
#'   expression profiles. Proceedings of the National Academy of Sciences,
#'   102(43), 15545-15550.
#'
#'   Liberzon, A., Subramanian, A., Pinchback, R., Thorvaldsdóttir, H., Tamayo,
#'   P., & Mesirov, J. P. (2011). Molecular signatures database (MSigDB) 3.0.
#'   Bioinformatics, 27(12), 1739-1740.
#'
#'   Liberzon, A., Birger, C., Thorvaldsdóttir, H., Ghandi, M., Mesirov, J. P.,
#'   & Tamayo, P. (2015). The molecular signatures database hallmark gene set
#'   collection. Cell systems, 1(6), 417-425.
#'
"hgsc"

#' An overlap network of signatures in the MSigDB
#'
#' The molecular signatures database (MSigDB) is a collection of over 25000 gene
#' expression signatures. Signatures in v7.2 are divided into 9 categories. This
#' network represents the precomputed overlap between all non-empty signatures
#' network represents the pre-computed overlap between all non-empty signatures
#' in the `msigdb` data. The network was computed using [computeMsigOverlap()]
#' and [computeMsigNetwork()] functions with default parameters.
#'
#' @format An igraph object, containing the overlap network computed using all
#'   non-empty gene expression signatures in MSigDB.
#' @docType data
#'
"msigOverlapNetwork"