test_that("word frequencies are computed correctly for edge cases", {
  #empty collections
  nullgsc = GSEABase::GeneSetCollection(list())
  freq = computeMsigWordFreq(nullgsc)

  expect_length(freq, 2)
  expect_equal(sapply(freq, ncol), c('Name' = 2, 'Short' = 2))
  expect_equal(sapply(freq, nrow), c('Name' = 0, 'Short' = 0))

  #non-empty gscs but empty gene sets
  emptygsc = GSEABase::GeneSetCollection(GSEABase::GeneSet())
  freq = computeMsigWordFreq(emptygsc)

  expect_length(freq, 2)
  expect_equal(sapply(freq, ncol), c('Name' = 2, 'Short' = 2))
  expect_equal(sapply(freq, nrow), c('Name' = 0, 'Short' = 0))

  emptygsc = GSEABase::GeneSetCollection(GSEABase::GeneSet(setName = 'a'),
                                         GSEABase::GeneSet(c('1', '2'), setName = 'b'))
  freq = computeMsigWordFreq(emptygsc)

  expect_length(freq, 2)
  expect_equal(sapply(freq, ncol), c('Name' = 2, 'Short' = 2))
  expect_equal(sapply(freq, nrow), c('Name' = 0, 'Short' = 0))
})

test_that("word frequencies (TF) are computed correctly", {
  #non-empty collections
  data(hgsc)
  estgsc = hgsc[grep('ESTROGEN', hgsc)]
  freq = computeMsigWordFreq(estgsc, measure = 'tf')

  expect_length(freq, 2)
  expect_equal(sapply(freq, nrow), c('Name' = 4, 'Short' = 5))
  expect_equal(sapply(freq, class), c('Name' = 'data.frame', 'Short' = 'data.frame'))
  expect_equal(freq$Name[freq$Name$word %in% c('early', 'late'), 2], c(1, 1))
  expect_equal(freq$Short[freq$Short$word %in% c('early', 'late'), 2], c(1, 1))
  expect_equal(freq$Name[freq$Name$word %in% c('estrogen'), 2], 2)
  expect_equal(freq$Short[freq$Short$word %in% c('estrogen'), 2], 2)
  expect_false('HALLMARK' %in% freq$Name$word)
  expect_false('to' %in% freq$Short$word)
})

test_that("word frequencies (TFIDF) are computed correctly", {
  #non-empty collections
  data(hgsc)
  estgsc = hgsc[grep('ESTROGEN', hgsc)]
  freq = computeMsigWordFreq(estgsc, measure = 'tfidf')

  expect_length(freq, 2)
  expect_equal(sapply(freq, nrow), c('Name' = 4, 'Short' = 5))
  expect_equal(sapply(freq, class), c('Name' = 'data.frame', 'Short' = 'data.frame'))
  expect_equal(freq$Name[freq$Name$word %in% c('early', 'late'), 2], c(1/3, 1/3))
  expect_equal(freq$Short[freq$Short$word %in% c('early', 'late'), 2], c(0.25, 0.25))
  expect_equal(freq$Name[freq$Name$word %in% c('estrogen'), 2], 0)
  expect_equal(freq$Short[freq$Short$word %in% c('estrogen'), 2], 0)
  expect_false('HALLMARK' %in% freq$Name$word)
  expect_false('to' %in% freq$Short$word)
})

test_that("blacklist works", {
  expect_true('blacklist' %in% getMsigBlacklist('blacklist'))
})
