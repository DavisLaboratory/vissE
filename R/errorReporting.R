checkGraph <- function(ig) {
  if (!is(ig, 'igraph')) {
    stop("'ig' is not an igraph object")
  }
}

checkNumericRange <- function(pvalue, pname, pmin = -Inf, pmax = Inf, eq = FALSE) {
  if (!is.numeric(pvalue) | length(pvalue) > 1)
    stop("'%s' should be a numeric of length 1", pvalue)

  if (eq) {
    inrange = pvalue >= pmin & pvalue <= pmax
    intstr = sprintf("[%s, %s]", pmin, pmax)
  } else {
    inrange = pvalue > pmin & pvalue < pmax
    intstr = sprintf("(%s, %s)", pmin, pmax)
  }

  if (!inrange) {
    if (is.infinite(pmin)) {
      stop(sprintf("'%s' should be < %s", pname, pmax))
    } else if (is.infinite(pmax)) {
      stop(sprintf("'%s' should be > %s", pname, pmin))
    } else {
      stop(sprintf("'%s' should be in the interval %s", pname, intstr))
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
  grpLen = sapply(groups, length)
  emptyGrps = names(groups)[grpLen == 0]
  if (length(emptyGrps) > 0) {
    emptyGrps = paste(emptyGrps, collapse = ', ')
    stop(sprintf("the following 'groups' contain no gene-sets: %s", emptyGrps))
  }

  #check for unknown gene-sets
  errorGrps = na.omit(sapply(names(groups), function(grpname) {
    if (!all(groups[[grpname]] %in% gscnames))
      return(grpname)
    return(NA)
  }))

  if (length(errorGrps) > 0) {
    errorGrps = paste(errorGrps, collapse = ', ')
    stop(sprintf("unknown gene-sets found in the following 'groups': %s", errorGrps))
  }
}

checkGenesetCollection <- function(gsc, pname) {
  if (!is(gsc, 'GeneSetCollection'))
    stop(sprintf("'%s' should be a GSEABase::GeneSetCollection object", pname))

  #check collection size
  if (length(gsc) == 0)
    stop(sprintf("'%s' cannot be an empty GeneSetCollection", pname))

  #check for empty gene-sets
  gscLen =  sapply(lapply(gsc, GSEABase::geneIds), length)
  emptyGscs = names(gsc)[gscLen == 0]
  if (length(emptyGscs) > 0) {
    emptyGscs = paste(emptyGscs, collapse = ', ')
    stop(sprintf("the following GeneSets in '%s' are empty: %s", pname, emptyGscs))
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
