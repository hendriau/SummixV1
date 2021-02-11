#' Example allele frequency data
#' 
#' reference data is 1000 Genomes and NAM. 
#' 1000 Genomes data was downloaded from ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/ on May 31, 2018
#' The IAM Affymetrix 6.0 data were downloaded from ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130711_native_american_admix_train (data accessed October 2018) and had been previously harmonized with the 1000 Genomes data.
#' Observed data is from gnomAD. gnomAD v2 data was downloaded from https://gnomad.broadinstitute.org/downloads on Oct. 11, 2018 
#' 
#' @docType data
#'
#' @usage data(ancestryData)
#'
#' @format Chromosome, SNP, base pair, reference and alternate alleles, 
#'  reference allele frequencies, observed allele frequencies
#'
#' @keywords datasets
#'
#' @examples
#' data("ancestryData")
"ancestryData"