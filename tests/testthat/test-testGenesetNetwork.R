library(GSEABase)
library(msigdb)

test_that("geneset intersection computation works", {
  g = paste0('ABCB', 1:12)
  gs1 = GeneSetCollection(
    GeneSet(g[1:5], setName = 'a', geneIdType = SymbolIdentifier()),
    GeneSet(g[4:8], setName = 'b', geneIdType = SymbolIdentifier())
  )
  gs2 = GeneSetCollection(
    GeneSet(g[6:9], setName = 'c', geneIdType = SymbolIdentifier()),
    GeneSet(g[9:12], setName = 'd', geneIdType = SymbolIdentifier())
  )
  id = msigdb::getMsigIdType(GeneSetCollection(gs1))

  expect_equal(intersectSize(geneIds(gs1), org = 'hs', idType = id)['a', 'b'], 2)
  expect_equal(intersectSize(geneIds(gs2), org = 'hs', idType = id)['c', 'd'], 1)
  expect_equal(intersectSize(geneIds(gs1), geneIds(gs2), org = 'hs', idType = id)['b', 'c'], 3)

  #works for single gene sets
  expect_equal(intersectSize(geneIds(gs1[2]), geneIds(gs2), org = 'hs', idType = id)['b', 'c'], 3)
})

test_that("geneset overlap computation works", {
  #empty collections
  nullgsc = GSEABase::GeneSetCollection(list())
  ov = computeMsigOverlap(nullgsc)
  expect_equal(nrow(ov), 0)
  expect_equal(ncol(ov), 3)
  expect_error(computeMsigOverlap(nullgsc, thresh = 2))
  expect_error(computeMsigOverlap(nullgsc, thresh = -2))
  expect_error(computeMsigOverlap(nullgsc, measure = 'euclidean'))

  #non-empty gscs but empty gene sets
  emptygsc = GSEABase::GeneSetCollection(GSEABase::GeneSet())
  expect_error(computeMsigOverlap(emptygsc, thresh = 1))

  emptygsc = GSEABase::GeneSetCollection(GSEABase::GeneSet(setName = 'a'),
                                         GSEABase::GeneSet(c('1', '2'), setName = 'b'))
  expect_error(computeMsigOverlap(emptygsc, thresh = 1))
  expect_error(computeMsigOverlap(emptygsc[2], emptygsc[1], thresh = 1))

  #2 geneset
  data("hgsc")
  estgsc = hgsc[grep('ESTROGEN', hgsc)]
  ov_jc = computeMsigOverlap(estgsc)
  ov_oc = computeMsigOverlap(estgsc, measure = 'ovlap')

  expect_equal(round(ov_jc[, 3], 3), 0.338)
  expect_equal(round(ov_oc[, 3], 3), 0.505)

  #hallmark geneset
  expect_equal(nrow(computeMsigOverlap(hgsc, nullgsc)), 0)
  expect_equal(ncol(computeMsigOverlap(hgsc, nullgsc)), 3)
  expect_equal(nrow(computeMsigOverlap(estgsc)), 1)
  expect_equal(nrow(computeMsigOverlap(hgsc, estgsc)), 4)
  expect_equal(nrow(computeMsigOverlap(hgsc)), 6)
  expect_equal(nrow(computeMsigOverlap(hgsc, thresh = 0)), 1225)
  expect_equal(nrow(computeMsigOverlap(hgsc, thresh = 1)), 0)
})

test_that("overlap network computation works", {
  #empty collections
  nullgsc = GSEABase::GeneSetCollection(list())
  ov = computeMsigOverlap(nullgsc)

  expect_error(computeMsigNetwork(ov, hgsc))

  #hallmark geneset
  data("hgsc")
  ov = computeMsigOverlap(hgsc)
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
  ov = computeMsigOverlap(hgsc)
  ig = computeMsigNetwork(ov, hgsc)
  grps = as.list(igraph::V(ig)$name)
  names(grps) = 1:length(grps)

  # expect_warning(plotMsigNetwork(ig, grps))
  expect_s3_class(plotMsigNetwork(ig, grps[1:3]), 'ggplot')
  expect_s3_class(plotMsigNetwork(ig, markGroups = list('a' = character())), 'ggplot')
  expect_error(plotMsigNetwork(ig, nodeSF = 0))
  expect_error(plotMsigNetwork(ig, nodeSF = -1))
  expect_error(plotMsigNetwork(ig, edgeSF = 0))
  expect_error(plotMsigNetwork(ig, edgeSF = -1))
  expect_error(plotMsigNetwork(ig, lytFunc = data.frame()))
  expect_error(plotMsigNetwork(ig, markGroups = character()))
  expect_error(plotMsigNetwork(ig, enrichStat = c(1)))
})

test_that("geneset neighbourhood works", {
  #hallmark geneset
  data("hgsc")
  ov = computeMsigOverlap(hgsc)
  ig = computeMsigNetwork(ov, hgsc)

  expect_equal(getMsigNeighbour('HALLMARK_ESTROGEN_RESPONSE_LATE', ig), 'HALLMARK_ESTROGEN_RESPONSE_EARLY')
  expect_length(getMsigNeighbour('HALLMARK_ESTROGEN_RESPONSE_LATE', ig, 0.5), 0)
  expect_error(getMsigNeighbour('FAKE_SIG', ig))
})
