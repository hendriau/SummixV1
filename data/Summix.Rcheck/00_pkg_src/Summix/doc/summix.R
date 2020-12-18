## ----eval=FALSE---------------------------------------------------------------
#  library("summix")
#  # load the data
#  data("ancestryData")
#  
#  # Estimate 5 reference ancestry proportion values for the gnomAD African/African Amercian ancestry group
#  summix( data = ancestryData,
#        reference=c("ref_AF_afr_1000G",
#            "ref_AF_eur_1000G",
#            "ref_AF_sas_1000G",
#            "ref_AF_iam_1000G",
#            "ref_AF_eas_1000G"),
#        observed="gnomad_AF_afr" )

## ----eval=FALSE---------------------------------------------------------------
#  library("summix")
#  data(ancestryData)
#  
#  head(ancestryData)
#  
#  tmp.aa<-adjAF(data   = ancestryData,
#       reference   = c("ref_AF_eur_1000G"),
#       observed    = "gnomad_AF_afr",
#       pi.target   = c(0, 1),
#       pi.observed = c(.15, .85))
#  tmp.aa$adjusted.AF[1:5,]

