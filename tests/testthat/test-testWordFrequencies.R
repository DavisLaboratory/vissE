test_that("word frequencies are computed correctly for edge cases", {
  #Null identifiers
  nullgsc = emptygsc = GSEABase::GeneSetCollection(GSEABase::GeneSet(setName = 'A'))
  expect_error(computeMsigWordFreq(nullgsc))
  
  #empty collections
  nullgsc = GSEABase::GeneSetCollection(list(), geneIdType = GSEABase::SymbolIdentifier())
  expect_error(computeMsigWordFreq(nullgsc))

  #non-empty gscs but empty gene sets
  emptygsc = GSEABase::GeneSetCollection(GSEABase::GeneSet(setName = 'A', geneIdType = GSEABase::SymbolIdentifier()))
  expect_error(computeMsigWordFreq(emptygsc))

  emptygsc = GSEABase::GeneSetCollection(GSEABase::GeneSet(setName = 'a', geneIdType = GSEABase::SymbolIdentifier()),
                                         GSEABase::GeneSet(c('ESR1', 'ESR2'), setName = 'b', geneIdType = GSEABase::SymbolIdentifier()))
  expect_error(computeMsigWordFreq(emptygsc))
})

test_that("word frequencies (TF) are computed correctly", {
  #non-empty collections
  data(hgsc)
  estgsc = hgsc[grep("ESTROGEN", hgsc)]
  freq = computeMsigWordFreq(estgsc, measure = 'tf')

  expect_length(freq, 2)
  expect_equal(sapply(freq, nrow), c('Name' = 4, 'Short' = 5))
  expect_equal(sapply(freq, class), c('Name' = 'data.frame', 'Short' = 'data.frame'))
  expect_equal(freq$Name[freq$Name$word %in% c('early', 'late'), 2], log(c(1, 1) + 1))
  expect_equal(freq$Short[freq$Short$word %in% c('early', 'late'), 2], log(c(1, 1) + 1))
  expect_equal(freq$Name[freq$Name$word %in% c('estrogen'), 2], log(2 + 1))
  expect_equal(freq$Short[freq$Short$word %in% c('estrogen'), 2], log(2 + 1))
  expect_false('HALLMARK' %in% freq$Name$word)
  expect_false('to' %in% freq$Short$word)
})

test_that("word frequencies (TFIDF) are computed correctly", {
  #non-empty collections
  data(hgsc)
  estgsc = hgsc[grep("ESTROGEN", hgsc)]
  idf = msigdb::getMsigdbIDF("hs", "7.2")
  freq = computeMsigWordFreq(estgsc, measure = 'tfidf', idf = idf)

  expect_length(freq, 2)
  expect_equal(sapply(freq, nrow), c('Name' = 4, 'Short' = 5))
  expect_equal(sapply(freq, class), c('Name' = 'data.frame', 'Short' = 'data.frame'))
  expect_equal(freq$Name[freq$Name$word %in% c("late"), 2], as.numeric(log(idf$Name["late"]) * log(1 + 1)), tolerance = 1e-5)
  expect_equal(freq$Short[freq$Short$word %in% c("early"), 2], as.numeric(log(idf$Short["early"]) * log(1 + 1)), tolerance = 1e-5)
  expect_equal(freq$Name[freq$Name$word %in% c("estrogen"), 2], as.numeric(log(idf$Name["estrogen"]) * log(2 + 1)), tolerance = 1e-5)
  expect_equal(freq$Short[freq$Short$word %in% c("estrogen"), 2], as.numeric(log(idf$Short["estrogen"]) * log(2 + 1)), tolerance = 1e-5)
  expect_false('HALLMARK' %in% freq$Name$word)
  expect_false('to' %in% freq$Short$word)
})

test_that("blacklist works", {
  expect_true('blacklist' %in% getMsigBlacklist('blacklist'))
})

test_that("weights work", {
  #non-empty collections
  data(hgsc)
  estgsc = hgsc[grep('ESTROGEN', hgsc)]
  
  w = c(2, 1)
  expect_error(computeMsigWordFreq(estgsc, w, measure = 'tf'))
  
  names(w) = names(estgsc)
  freq = computeMsigWordFreq(estgsc, w, measure = 'tf')
  expect_equal(freq$Name[freq$Name$word %in% c('early', 'late'), 2], as.numeric(log(c(1, 1) * w + 1)))
  expect_equal(freq$Short[freq$Short$word %in% c('early', 'late'), 2], as.numeric(log(c(1, 1) * w + 1)))
})
