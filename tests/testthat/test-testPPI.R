test_that("PPI network gets computed", {
  data('hgsc')
  ppihs = getIMEX('hs')
  ppimm = getIMEX('mm')
  grps = list('early' = 'HALLMARK_ESTROGEN_RESPONSE_EARLY', 'late' = 'HALLMARK_ESTROGEN_RESPONSE_LATE')
  
  expect_s3_class(plotMsigPPI(ppihs, hgsc, grps), 'ggplot')
  expect_error(plotMsigPPI(ppimm, hgsc, grps), 'organism')
})
