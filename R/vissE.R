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

#' Binary membership matrix for the Human MSigDB
#'
#' This object stores the Human molecular signatures database (MSigDB) in binary
#' format as a membership matrix. Gene signatures are along the rows and Entrez
#' IDs are along the columns.
#' @format A dgCMatrix (sparse) object, with gene sets along the rows and Entrez
#'   IDs along the columns.
#' @docType data
#'   
"mem_mat_hs"

#' Binary membership matrix for the Mouse MSigDB
#'
#' This object stores the Mouse molecular signatures database (MSigDB) in binary
#' format as a membership matrix. Gene signatures are along the rows and Entrez
#' IDs are along the columns.
#' @format A dgCMatrix (sparse) object, with gene sets along the rows and Entrez
#'   IDs along the columns.
#' @docType data
#'   
"mem_mat_mm"

.myDataEnv <- new.env(parent = emptyenv()) # not exported

.data_internal <- function(dataset) {
  if (!exists(dataset, envir = .myDataEnv)) {
    utils::data(list = c(dataset), envir = .myDataEnv)
  }
}
