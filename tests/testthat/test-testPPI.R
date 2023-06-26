test_that("PPI network gets computed", {
  library(msigdb)
  data('hgsc')
  ppihs = getIMEX("hs", version = "2022-02-03")
  ppimm = getIMEX("mm", version = "2022-02-03")
  grps = list('early' = 'HALLMARK_ESTROGEN_RESPONSE_EARLY', 'late' = 'HALLMARK_ESTROGEN_RESPONSE_LATE')
  
  expect_s3_class(plotMsigPPI(ppihs, hgsc, grps), 'ggplot')
  expect_error(plotMsigPPI(ppimm, hgsc, grps), 'organism')
})
