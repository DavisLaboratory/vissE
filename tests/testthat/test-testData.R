test_that("hgsc data is in the correct format", {
  library(GSEABase)

  data(hgsc)

  expect_true(all(sapply(lapply(hgsc, geneIds), length) > 0))
  expect_true(all(sapply(lapply(hgsc, collectionType), class) %in% 'BroadCollection'))
  expect_true(all(sapply(lapply(hgsc, collectionType), bcCategory) %in% 'h'))
  expect_length(hgsc, 50)
})
