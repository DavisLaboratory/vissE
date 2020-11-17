computeMsigWordFreq <- function(msigGsc, rmwords = getMsigBlacklist()) {
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
  # Text stemming
  # docs <- tm_map(docs, stemDocument)

  #compute frequencies
  dtms = lapply(docs, tm::TermDocumentMatrix)
  v = lapply(dtms, function(x) sort(apply(x, 1, sum), decreasing = TRUE))
  d = lapply(v, function(x) data.frame(word = names(x), freq = x))

  #remove GSE labels
  d = lapply(d, function(x) x[!grepl('gse', x$word), ])

  if (is.null(nrow(d$Name))) {
    d$Name = data.frame('word' = character(), 'freq' = numeric())
  }
  if (is.null(nrow(d$Short))) {
    d$Short = data.frame('word' = character(), 'freq' = numeric())
  }

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
#'
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
    'reactome',
    'up'
  )

  #add custom words
  rmwords = sort(unique(c(rmwords, custom)))

  return(rmwords)
}
