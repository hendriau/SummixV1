test_that("adjAF() works on ancestryData", {
  expect <- c(0.54425727, 0.05379769, 0.13027010, 0.93017307, 0.81593620)
  data("ancestryData")
  testval <- adjAF(data   = ancestryData,
                   reference   = c("ref_AF_eur_1000G"),
                   observed    = "gnomad_AF_afr",
                   pi.target   = c(0, 1),
                   pi.observed = c(.15, .85))
  expect_equal(testval$adjusted.AF[1:5,"adjustedAF"], expected = expect)
})
