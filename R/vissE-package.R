#' @details
#' This package supports four workflows to enhance gene set enrichment analysis:
#' \enumerate{
#'   \item Clustering results from a gene set enrichment analysis (e.g. using
#'   limma::fry, singscore or GSEA). The functions required for this analysis
#'   are \code{\link{computeMsigOverlap}}, \code{\link{computeMsigNetwork}} and
#'   \code{\link{plotMsigNetwork}}.
#'   \item Interpreting gene set clusters (identified in the first analysis) by
#'   performing text-mining of gene set names and descriptions. The main
#'   function required to perform text-mining of gene sets is
#'   \code{\link{plotMsigWordcloud}}. Other functions can be used to access
#'   intermmediate results.
#'   \item Visualise gene-level statistics for gene set clusters identified in
#'   the first analysis to link back gene set clusters to the genes of interest.
#'   This can be done using the \code{\link{plotGeneStats}} function.
#'   \item Identifying gene sets similar to a list of genes identified from a DE
#'   analysis using set overlap measures. This can be done using the
#'   \code{\link{characteriseGeneset}} function.
#' }
#' 
#' @keywords internal
"_PACKAGE"

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL