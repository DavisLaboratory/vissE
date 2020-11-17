test_that("geneset overlap computation works", {
  #empty collections
  nullgsc = GSEABase::GeneSetCollection(list())
  ov = computeMsigOverlap(nullgsc)
  expect_equal(nrow(ov), 0)
  expect_equal(ncol(ov), 3)
  expect_error(computeMsigOverlap(nullgsc, thresh = 2))
  expect_error(computeMsigOverlap(nullgsc, thresh = -2))
  expect_error(computeMsigOverlap(nullgsc, measure = 'euclidean'))

  #2 geneset
  data(hgsc)
  estgsc = hgsc[grep('ESTROGEN', hgsc)]
  ov_jc = computeMsigOverlap(estgsc)
  ov_oc = computeMsigOverlap(estgsc, measure = 'ovlap')

  expect_equal(round(ov_jc[, 3], 3), 0.338)
  expect_equal(round(ov_oc[, 3], 3), 0.505)

  #hallmark geneset
  expect_equal(nrow(computeMsigOverlap(hgsc)), 14)
  expect_equal(nrow(computeMsigOverlap(hgsc, thresh = 0)), 1225)
  expect_equal(nrow(computeMsigOverlap(hgsc, thresh = 1)), 0)
})

test_that("overlap network computation works", {
  #empty collections
  nullgsc = GSEABase::GeneSetCollection(list())
  ov = computeMsigOverlap(nullgsc)

  expect_error(computeMsigNetwork(ov, hgsc))

  #hallmark geneset
  data(hgsc)
  ov = computeMsigOverlap(hgsc)
  ig = computeMsigNetwork(ov, hgsc)

  expect_equal(class(ig), 'igraph')
  expect_length(igraph::V(ig), 21)
  expect_length(igraph::E(ig), 14)
  expect_error(computeMsigNetwork(ov, hgsc[1:10]))
  expect_error(computeMsigNetwork(rbind(ov, c('a', 'a','1')), hgsc[1:10]))
})