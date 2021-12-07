library(GSEABase)
library(msigdb)

test_that("geneset overlap computation works", {
  #empty collections
  nullgsc = GeneSetCollection(list())
  expect_error(computeMsigOverlap(nullgsc, thresh = 0.15, measure = 'jaccard'))
  
  #non-empty gscs but empty gene sets
  emptygsc = GeneSetCollection(GeneSet(setName = 'A', geneIdType = SymbolIdentifier()))
  expect_error(computeMsigOverlap(emptygsc, thresh = 1, measure = 'jaccard'))
  
  emptygsc = GSEABase::GeneSetCollection(GSEABase::GeneSet(setName = 'a'),
                                         GSEABase::GeneSet(c('1', '2'), setName = 'b'))
  expect_error(computeMsigOverlap(emptygsc, thresh = 1, measure = 'jaccard'))
  expect_error(computeMsigOverlap(emptygsc[2], emptygsc[1], thresh = 1, measure = 'jaccard'))
  
  #2 geneset
  data("hgsc")
  estgsc = hgsc[grep('ESTROGEN', hgsc)]
  ov_ari = computeMsigOverlap(estgsc, thresh = -Inf, measure = 'ari')
  ov_jc = computeMsigOverlap(estgsc, thresh = 0.15, measure = 'jaccard')
  ov_oc = computeMsigOverlap(estgsc, thresh = 0.15, measure = 'ovlap')
  
  expect_equal(round(ov_ari[, 3], 3), 0.091)
  expect_equal(round(ov_jc[, 3], 3), 0.338)
  expect_equal(round(ov_oc[, 3], 3), 0.505)
  
  #hallmark geneset
  expect_error(nrow(computeMsigOverlap(hgsc, nullgsc, thresh = 0.15, measure = 'jaccard')))
  expect_error(ncol(computeMsigOverlap(hgsc, nullgsc, thresh = 0.15, measure = 'jaccard')))
  expect_equal(nrow(computeMsigOverlap(estgsc, thresh = 0.15, measure = 'jaccard')), 1)
  expect_equal(nrow(computeMsigOverlap(hgsc, estgsc, thresh = 0.15, measure = 'jaccard')), 4)
  expect_equal(nrow(computeMsigOverlap(hgsc, thresh = 0.15, measure = 'jaccard')), 6)
  expect_equal(nrow(computeMsigOverlap(hgsc, thresh = 0, measure = 'jaccard')), 1225)
  expect_equal(nrow(computeMsigOverlap(hgsc, thresh = 1, measure = 'jaccard')), 0)
})

test_that("overlap network computation works", {
  #hallmark geneset
  data("hgsc")
  ov = computeMsigOverlap(hgsc, thresh = 0.15, measure = 'jaccard')
  ig = computeMsigNetwork(ov, hgsc)
  
  expect_s3_class(ig, 'igraph')
  expect_length(igraph::V(ig), 12)
  expect_length(igraph::E(ig), 6)
  expect_error(computeMsigNetwork(ov, hgsc[1:10]))
  expect_error(computeMsigNetwork(rbind(ov, c('a', 'a','1')), hgsc[1:10]))
})

test_that("plot functions work", {
  #hallmark geneset
  data("hgsc")
  ov = computeMsigOverlap(hgsc, thresh = 0.15, measure = 'jaccard')
  ig = computeMsigNetwork(ov, hgsc)
  grps = as.list(igraph::V(ig)$name)
  names(grps) = 1:length(grps)
  
  # expect_warning(plotMsigNetwork(ig, grps))
  expect_s3_class(plotMsigNetwork(ig, grps[1:3]), 'ggplot')
  expect_error(plotMsigNetwork(ig, markGroups = list('a' = character())), 'groups')
  expect_error(plotMsigNetwork(ig, nodeSF = 0))
  expect_error(plotMsigNetwork(ig, nodeSF = -1))
  expect_error(plotMsigNetwork(ig, edgeSF = 0))
  expect_error(plotMsigNetwork(ig, edgeSF = -1))
  expect_error(plotMsigNetwork(ig, lytFunc = data.frame()))
  expect_error(plotMsigNetwork(ig, markGroups = character()))
  expect_error(plotMsigNetwork(ig, genesetStat = c(1)))
})

test_that("geneset neighbourhood works", {
  #hallmark geneset
  data("hgsc")
  ov = computeMsigOverlap(hgsc, thresh = 0.15, measure = 'jaccard')
  ig = computeMsigNetwork(ov, hgsc)
  
  expect_equal(getMsigNeighbour('HALLMARK_ESTROGEN_RESPONSE_LATE', ig), 'HALLMARK_ESTROGEN_RESPONSE_EARLY')
  expect_length(getMsigNeighbour('HALLMARK_ESTROGEN_RESPONSE_LATE', ig, 0.5), 0)
  expect_error(getMsigNeighbour('FAKE_SIG', ig))
})
