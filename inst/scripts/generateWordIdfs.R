library(org.Hs.eg.db)
library(org.Mm.eg.db)

computeIdf <- function(msigGsc) {
  rmwords = vissE:::getMsigBlacklist()
  
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
  
  #compute idf
  dtms = lapply(docs, tm::TermDocumentMatrix)
  dtms = lapply(dtms, as.matrix)
  
  #compute IDF
  idfs = lapply(dtms, function(x) {
    idf = log(ncol(x) / rowSums(x != 0))
    return(idf)
  })
  #sort names to quicken searches
  idfs = lapply(idfs, function(x) {
    x[sort(names(x))]
  })
  
  return(idfs)
}

msigdb.hs = msigdb::msigdb.v7.2.hs.SYM()
msigdb.mm = msigdb::msigdb.v7.2.mm.SYM()

msigdb.hs = msigdb::appendKEGG(msigdb.hs)
msigdb.mm = msigdb::appendKEGG(msigdb.mm)

idf_hs = computeIdf(msigdb.hs)
idf_mm = computeIdf(msigdb.mm)

#----namemap for membership matrix----
selcolc = c('h', 'c5')
mem_mat_hs_map = unique(unlist(geneIds(subsetCollection(msigdb.hs, selcolc))))
mem_mat_mm_map = unique(unlist(geneIds(subsetCollection(msigdb.mm, selcolc))))

usethis::use_data(
  idf_hs,
  idf_mm,
  mem_mat_hs_map,
  mem_mat_mm_map,
  internal = TRUE,
  overwrite = TRUE
)

