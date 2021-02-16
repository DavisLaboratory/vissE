#' Compute word frequencies for a single MSigDB collection
#'
#' @param msigGsc a GeneSetCollection object, containing gene sets from the
#'   MSigDB. The [GSEABase::getBroadSets()] function can be used to parse XML
#'   files downloaded from MSigDB.
#' @param measure a character, specifying how frequencies should be computed.
#'   "tf" uses term frequencies and "tfidf" (default) applies inverse document
#'   frequency weights to term frequencies.
#' @param rmwords a character vector, containing a blacklist of words to discard
#'   from the analysis.
#'
#' @return a list, containing two data.frames summarising the results of the
#'   frequency analysis on gene set names and short descriptions.
#' @export
#'
#' @examples
#' data(hgsc)
#' freq <- computeMsigWordFreq(hgsc, measure = 'tfidf')
#'
computeMsigWordFreq <- function(msigGsc, measure = c('tfidf', 'tf'), rmwords = getMsigBlacklist()) {
  measure = match.arg(measure)
  stopifnot('GeneSetCollection' %in% class(msigGsc))

  signames = sapply(msigGsc, GSEABase::setName)
  sigdesc_s = sapply(msigGsc, GSEABase::description)
  docs = list('Name' = signames, 'Short' = sigdesc_s)
  docs = lapply(docs, unique)

  #text-mining
  docs = lapply(docs, function(d) tm::Corpus(tm::VectorSource(d)))
  toSpace <- tm::content_transformer(function (x, pattern) gsub(pattern, " ", x))
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

  #compute frequencies
  dtms = lapply(docs, tm::TermDocumentMatrix)
  dtms = lapply(dtms, tm::weightTf)
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
  org = msigdb::getMsigOrganism(msigGsc, idType)
  if (org %in% 'hs') {
    idf = idf_hs
  } else {
    idf = idf_mm
  }
  
  #compute tfidfs
  d = lapply(d, function(x) {
    x$freq = log(1 + x$freq)
    return(x)
  })
  
  if (measure %in% 'tfidf'){
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

#' Blacklist words for MSigDB gene set text mining
#'
#' List of words to discard when performing text mining MSigDB gene set names
#' and short descriptions.
#'
#' @param custom a character vector, containing list of words to add onto
#'   existing blacklist.
#'
#' @return a character vector, containing list of blacklist works
#' @export
#'
#' @examples
#' getMsigBlacklist('blacklist')
#'
getMsigBlacklist <- function(custom = c()) {
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
    'gse',
    'hallmark',
    'kegg',
    'module',
    'morf',
    'neighborhood',
    'pathway',
    'reactome',
    'up'
  )

  #add custom words
  rmwords = sort(unique(c(rmwords, custom)))

  return(rmwords)
}
