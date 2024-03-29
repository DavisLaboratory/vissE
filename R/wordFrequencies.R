#' Compute word frequencies for a single MSigDB collection
#'
#' @param msigGsc a GeneSetCollection object, containing gene sets from the
#'   MSigDB. The [GSEABase::getBroadSets()] function can be used to parse XML
#'   files downloaded from MSigDB.
#' @param measure a character, specifying how frequencies should be computed.
#'   "tf" uses term frequencies and "tfidf" (default) applies inverse document
#'   frequency weights to term frequencies.
#' @param rmwords a character vector, containing an exclusion list of words to discard
#'   from the analysis.
#' @param weight a named numeric vector, containing weights to apply to each
#'   gene-set. This can be -log10(FDR), -log10(p-value) or an enrichment score
#'   (ideally unsigned).
#' @param version a character, specifying the version of msigdb to use (see
#'   `msigdb::getMsigdbVersions()`).
#' @param org a character, specifying the organism to use. This can either be
#'   "auto" (default), "hs" or "mm".
#' @param idf a list of named numeric vectors, specifying inverse document frequencies to use to penalise terms from gene-set names and short descriptions. This should be a vector of length 2 with names "Name" and "Short". Numeric vectors should contain weights and names should represent the term. Precomputed versions can be retrieved using the [msigdb::getMsigdbIDF()].
#'
#' @return a list, containing two data.frames summarising the results of the
#'   frequency analysis on gene set names and short descriptions.
#' @export
#'
#' @examples
#' data(hgsc)
#' freq <- computeMsigWordFreq(hgsc, measure = 'tfidf')
#' 
computeMsigWordFreq <-
  function(msigGsc,
           weight = NULL,
           measure = c('tfidf', 'tf'),
           version = msigdb::getMsigdbVersions(),
           org = c('auto', 'hs', 'mm'),
           rmwords = getMsigExclusionList(),
           idf = NULL) {
  #check params  
  measure = match.arg(measure)
  version = match.arg(version)
  org = match.arg(org)
  checkGenesetCollection(msigGsc, 'msigGsc')
  if (!is.null(weight)) {
    if (any(weight <= 0))
      stop("all weights in 'weight' should be greater than 0")
    if (is.null(names(weight)))
      stop("'weight' should be a named vector")
    missingWeights = setdiff(names(msigGsc), names(weight))
    if (length(missingWeights) > 0) {
      missingWeights = paste(missingWeights, collapse = ', ')
      stop('the following gene-sets are missing weights: %s', missingWeights)
    }
  }

  #check weights
  if (is.null(weight)) {
    weight = rep(1, length(msigGsc))
    names(weight) = sapply(msigGsc, GSEABase::setName)
  } else{
    weight = weight[names(msigGsc)]
  }
  
  #extract text data from signatures
  signames = sapply(msigGsc, GSEABase::setName)
  sigdesc_s = sapply(msigGsc, GSEABase::description)
  docs = list('Name' = signames, 'Short' = sigdesc_s)
  docs = lapply(docs, unique)

  #text-mining
  docs = lapply(docs, function(d) tm::Corpus(tm::VectorSource(d)))
  toSpace <- tm::content_transformer(function (x, pattern) gsub(pattern, " ", x))
  suppressWarnings({
    docs = lapply(docs, function(d) tm::tm_map(d, toSpace, "_"))
    docs = lapply(docs, function(d) tm::tm_map(d, toSpace, "/"))
    docs = lapply(docs, function(d) tm::tm_map(d, toSpace, "@"))
    docs = lapply(docs, function(d) tm::tm_map(d, toSpace, "\\|"))
    docs = lapply(docs, function(d) tm::tm_map(d, toSpace, "\\("))
    docs = lapply(docs, function(d) tm::tm_map(d, toSpace, "\\)"))
    
    # Convert the text to lower case
    docs = lapply(docs, function(d) tm::tm_map(d, tm::content_transformer(tolower)))
    # Remove numbers
    # docs = lapply(docs, function(d) tm_map(d, removeNumbers))
    # Remove english common stopwords
    docs = lapply(docs, function(d) tm::tm_map(d, tm::removeWords, tm::stopwords('english')))
    # Remove your own stop word
    # specify your stopwords as a character vector
    docs = lapply(docs, function(d) tm::tm_map(d, tm::removeWords, rmwords))
    # Remove punctuations
    docs = lapply(docs, function(d) tm::tm_map(d, tm::removePunctuation))
    # Eliminate extra white spaces
    docs = lapply(docs, function(d) tm::tm_map(d, tm::stripWhitespace))
    # Remove full numbers
    docs = lapply(docs, function(d) tm::tm_filter(d, function(x) !grepl('\\b[0-9]+\\b', x)))
    # Text lemmatisation
    docs = lapply(docs, function(d) tm::tm_map(d, textstem::lemmatize_strings))
  })

  #compute frequencies
  dtms = lapply(docs, tm::TermDocumentMatrix)
  dtms = lapply(dtms, tm::weightTf)
  #apply weights
  dtms = lapply(dtms, function(x) {
    #log for TF-IDF computation
    x$v = x$v * weight[x$j]
    return(x)
  })
  v = lapply(dtms, function(x) apply(x, 1, sum))
  d = lapply(v, function(x) data.frame(word = names(x), freq = x))

  #remove GSE labels
  d = lapply(d, function(x) x[!grepl('gse', x$word), ])

  if (is.null(nrow(d$Name))) {
    d$Name = data.frame('word' = character(), 'freq' = numeric())
  }
  if (is.null(nrow(d$Short))) {
    d$Short = data.frame('word' = character(), 'freq' = numeric())
  }
  if (length(msigGsc) == 0)
    return(d)
  
  #identify the organism and use the correct IDFs
  idType = msigdb::getMsigIdType(msigGsc)
  if (org %in% 'auto') {
    org = msigdb::getMsigOrganism(msigGsc, idType)
  }
  
  #compute log TF
  d = lapply(d, function(x) {
    x$freq = log(1 + x$freq)
    return(x)
  })
  
  #compute tfidfs
  if (measure %in% 'tfidf'){
    #select pre-processed data for organism
    if (is.null(idf)) {
      idf = msigdb::getMsigdbIDF(org, version)
    } else {
      if (
        !all(c('Name', 'Short') %in% names(idf)) |
        !all(sapply(idf, \(x) !is.null(names(x)))) |
        !all(sapply(idf, is.numeric))
      ) {
        stop('Invalid precomputed IDF object provided')
      }
    }

    d = mapply(function(x, i) {
      x$freq = log(i[x$word]) * x$freq
      return(x)
    }, d, idf, SIMPLIFY = FALSE)
  }
  
  #sort by TF-IDF
  d = lapply(d, function(x) {
    x = x[order(x$freq, decreasing = TRUE), ]
    x = x[!is.na(x$freq), ]
    rownames(x) = NULL
    return(x)
  })

  return(d)
}

#' Exclusion list of words for MSigDB gene set text mining
#'
#' List of words to discard when performing text mining MSigDB gene set names
#' and short descriptions.
#'
#' @param custom a character vector, containing list of words to add onto
#'   existing exclusion list.
#'
#' @return a character vector, containing words to be excluded from the text mining analysis.
#' @export
#'
#' @examples
#' getMsigExclusionList('remove')
#'
getMsigExclusionList <- function(custom = c()) {
  rmwords = c(
    'biocarta',
    'car',
    'cells',
    'dn',
    'expression',
    'gcm',
    'gene',
    'genes',
    'gnf2',
    'gobp',
    'gomf',
    'gocc',
    'gse',
    'hallmark',
    'kegg',
    'module',
    'morf',
    'neighborhood',
    'pathway',
    'pid',
    'reactome',
    'up',
    'wp'
  )

  #add custom words
  rmwords = sort(unique(c(rmwords, custom)))

  return(rmwords)
}
