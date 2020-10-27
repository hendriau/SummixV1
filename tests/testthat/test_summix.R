test_that("summix() works on ancestryData", {
  data("ancestryData")
  testval <- summix( data = ancestryData, 
                     reference=c("ref_AF_afr_1000G", 
                                 "ref_AF_eur_1000G", 
                                 "ref_AF_sas_1000G", 
                                 "ref_AF_iam_1000G", 
                                 "ref_AF_eas_1000G"), 
                     observed="gnomad_AF_afr", 
                     pi.start = c(.2, .2, .2, .2, .2) )
  
  expect_equal(length(testval), expected = 9)
})
