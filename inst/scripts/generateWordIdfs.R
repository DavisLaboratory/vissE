library(org.Hs.eg.db)
library(org.Mm.eg.db)

computeIdf <- function(msigGsc) {
  # msigGsc = msigdb::msigdb.hs.SYM()
  msigGsc = msigdb::appendKEGG(msigGsc)
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
  # Text stemming
  dicts = docs
  docs = lapply(docs, function(d) tm::tm_map(d, tm::stemDocument, language = 'english'))
  
  #compute idf
  dtms = lapply(docs, tm::TermDocumentMatrix)
  dtms = lapply(dtms, as.matrix)
  stemmap = mapply(function(x, d) {
    smap = tm::stemCompletion(rownames(x), d, type = 'prevalent')
    names(smap) = rownames(x)
    return(smap)
  }, dtms, dicts, SIMPLIFY = FALSE)
  
  #combine frequencies for completed-stemmed words
  dtms = mapply(function(x, smap) {
    rowsum(x, as.factor(smap[rownames(x)]))
  }, dtms, stemmap)
  
  idfs = lapply(dtms, function(x) {
    idf = log(ncol(x) / rowSums(x != 0))
    return(idf)
  })
  
  #create a df
  df_name = merge(
    data.frame('Stem' = names(stemmap$Name), 'Complete' = stemmap$Name),
    data.frame('Complete' = names(idfs$Name), 'Name.IDF' = idfs$Name),
    all = TRUE
  )
  df_short = merge(
    data.frame('Stem' = names(stemmap$Short), 'Complete' = stemmap$Short),
    data.frame('Complete' = names(idfs$Short), 'Short.IDF' = idfs$Short),
    all = TRUE
  )
  df_idf = merge(df_name, df_short, all = TRUE)
  df_idf = df_idf[df_idf$Complete != '', ]
  rownames(df_idf) = NULL
  
  return(df_idf)
}

e = new.env()
load('../msigdb/msigdb.hs.SYM.rda', envir = e)
load('../msigdb/msigdb.mm.SYM.rda', envir = e)

df_idf_hs = computeIdf(e$msigdb.hs.SYM)
df_idf_mm = computeIdf(e$msigdb.mm.SYM)

usethis::use_data(df_idf_hs, internal = TRUE)
usethis::use_data(df_idf_mm, internal = TRUE)

