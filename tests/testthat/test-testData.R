test_that("hgsc data is in the correct format", {
  library(GSEABase)

  data(hgsc)

  expect_true(all(sapply(lapply(hgsc, geneIds), length) > 0))
  expect_true(all(sapply(lapply(hgsc, collectionType), class) %in% 'BroadCollection'))
  expect_true(all(sapply(lapply(hgsc, collectionType), bcCategory) %in% 'h'))
  expect_length(hgsc, 50)
})

test_that("msigOverlapNetwork data is in the correct format", {
  # library(GSEABase)
  # library(igraph)
  # library(msigdb)
  # 
  # msigdb = msigdb.hs.SYM()
  # data(msigOverlapNetwork)
  # 
  # expect_true(all(V(msigOverlapNetwork)$name %in% sapply(msigdb, setName)))
})
