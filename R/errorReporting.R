checkGraph <- function(ig) {
  if (!is(ig, 'igraph')) {
    stop("'ig' is not an igraph object")
  }
}

checkNumericRange <- function(pvalue, pname, pmin = -Inf, pmax = Inf) {
  if (pvalue < pmin | pvalue > pmax) {
    if (is.infinite(pmin)) {
      stop(sprintf("'%s' should be < %s", pname, pmax))
    } else if (is.infinite(pmax)) {
      stop(sprintf("'%s' should be > %s", pname, pmin))
    } else {
      stop(sprintf("'%s' should be in the interval (%s, %s)", pname, pmin, pmax))
    }
  }
}

checkGroups <- function(groups, gscnames) {
  if (!is.list(groups))
    stop("'groups' should be a list")
  
  if (length(groups) == 0)
    stop("'groups' should be a non-empty list")
  
  if (is.null(names(groups)))
    stop("'groups' must be a named list")
  
  #check for empty groups
  grp_len = sapply(groups, length)
  empty_grps = names(groups)[grp_len == 0]
  if (length(empty_grps) > 0) {
    empty_grps = paste(empty_grps, collapse = ', ')
    stop(sprintf("the following 'groups' contain no gene-sets: %s", empty_grps))
  }
  
  #check for unknown gene-sets
  lapply(names(groups), function(grpname) {
    if (!all(groups[[grpname]] %in% gscnames))
      stop(sprintf("unknown gene-sets found in the following 'groups':", grpname))
  })
}

checkGenesetCollection <- function(gsc, pname) {
  if (!is(gsc, 'GeneSetCollection'))
    stop(sprintf("'%s' should be a GeneSetCollection object", pname))
  
  #check collection size
  if (length(gsc) == 0)
    stop(sprintf("'%s' cannot be an empty GeneSetCollection", pname))
  
  #check for empty gene-sets
  grp_len =  sapply(lapply(gsc, GSEABase::geneIds), length)
  empty_gscs = names(gsc)[grp_len == 0]
  if (length(empty_gscs) > 0) {
    empty_gscs = paste(empty_gscs, collapse = ', ')
    stop(sprintf("the following GeneSets in '%s' are empty: %s", pname, empty_grps))
  }
}

checkGeneStat <- function(geneStat) {
  if (is.null(names(geneStat)))
    stop("'geneStat' should be a named vector")
  
  if (!is.numeric(geneStat))
    stop("'geneStat' should be a numeric vector")
}

checkGenesetStat <- function(genesetStat) {
  if (is.null(names(genesetStat)))
    stop("'genesetStat' should be a named vector")
  
  if (!is.numeric(genesetStat))
    stop("'genesetStat' should be a numeric vector")
}
