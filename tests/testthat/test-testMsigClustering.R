library(GSEABase)
library(msigdb)

test_that("geneset overlap computation works", {
  #hallmark geneset
  data("hgsc")
  ov_ari = computeMsigOverlap(hgsc, thresh = 0.25, measure = 'ari')
  ov_ig = computeMsigNetwork(ov_ari, hgsc)
  
  grps = findMsigClusters(ov_ig)
  expect_true(any(sapply(grps, function(x) all(grepl('ESTROGEN', x)))))
})