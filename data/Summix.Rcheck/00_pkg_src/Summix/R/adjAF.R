#' Ancestry adjusted allele frequencies
#' 
#' Adjusts allele frequencies for heterogeneous populations in genetic data given proportion of reference ancestry groups
#'
#' @param data dataframe of unadjusted allele frequency for observed group, K-1 reference ancestry allele frequencies for N SNPs
#' @param reference character vector of the column names for K-1 reference ancestry groups. The name of the last reference ancestry group is not included as that group is not used to estimate the adjusted allele frequencies. 
#' @param observed character value for the column name of observed data group
#' @param pi.target numeric vector of the mixture proportions for K reference ancestry groups in the target sample or subject. The order must match the order of the reference columns with the last entry matching the missing reference group.
#' @param pi.observed  numeric vector of the mixture proportions for K reference ancestry groups for the observed group. The order must match the order of the reference columns with the last entry matching the missing reference group.
#' @return pi: table of input reference ancestry groups, pi.observed, and pi.target
#' @return observed.data: name of the data column for the observed group from which adjusted ancestry allele frequency is estimated
#' @return Nsnps: number of SNPs for which adjusted AF is estimated
#' @return adjusted.AF: data frame of original data with an appended column of adjusted allele frequencies
#'
#' @author Gregory Matesi, \email{gregory.matesi@ucdenver.edu}
#' @author Audrey Hendricks, \email{audrey.hendricks@ucdenver.edu}
#'
#' @reference
#' @keywords genetics mixture ancestry
#' 
#' @seealso \code{\link{summix}} for estimating the proportion of reference ancestry groups and \url{https://github.com/hendriau/Summix} for further documentation
#' 
#' @examples
#' data(ancestryData)
#' tmp.aa<-adjAF(data   = ancestryData,
#'     reference   = c("ref_AF_eur_1000G"),
#'     observed    = "gnomad_AF_afr",
#'     pi.target   = c(0, 1),
#'     pi.observed = c(.15, .85))
#' tmp.aa$adjusted.AF[1:5,]
#' @export


# INPUTS:
# data (dataframe):
#     CHR  RSID       POS      REF ALT  ref_eur     ...  ref_iam   target_afr obs_amr    obs_oth
#     1    rs2887286  1156131  C  T   0.173275495 ...  0.7093    0.4886100  0.52594300 0.22970500
#     1    rs41477744 2329564  A  G   0.001237745 ...  0.0000    0.0459137  0.00117925 0.00827206
#     1    rs9661525  2952840  G  T   0.168316089 ...  0.2442    0.1359770  0.28605200 0.15561700
#     1    rs2817174  3044181  C  T   0.428212624 ...  0.5000    0.8548790  0.48818000 0.47042500
#     1    rs12139206 3504073  T  C   0.204214851 ...  0.3372    0.7241780  0.29550800 0.25874800
#     1    rs7514979  3654595  T  C   0.004950604 ...  0.0000    0.3362490  0.01650940 0.02481620
#
# reference: Names of K-1 columns with reference AF. The name of the last reference group is not included as the AFs are not used for this group. Ex. "gnomad_AF_nfe" (note, African ancestry is excluded)
#
# observed: name of column with observed AF to update to ancestry adjusted. Ex. "gnomad_AF_afr"
#
# pi.target: proportion of ancestry for K ancestry groups in the target sample. The order must match the order of the reference columns with the last entry matching the missing group in reference. Ex. c(0,1) for 0 European ancestry and 1 African ancestry 
#
# pi.observed: Proportion of ancestry for K ancestry groups in the observed data (must match the observed column). The order must match the order of the reference columns with the last entry matching the missing group in reference. Ex. c(0.16,0.84) for 0.16 European ancestry and 0.84 African ancestry 
#
#
# 
#   
#

adjAF <- function(
                      data         ,
                      reference    ,
                      observed     ,
                      pi.target    , 
                      pi.observed  ){
  
  
  if(length(pi.target) != length(pi.observed)){
    stop("ERROR Please make sure that the lengths of pi.target and pi.observed match.")
  }
  if(length(reference) != length(pi.target)-1){
    stop("ERROR: Please make sure that you are only entering k-1 ancestries in reference.")
  }

  
  if(!is(object = data, class2 = "data.frame")){
    stop("ERROR: data must be a data.frame as described in the vignette")
  }

  if(typeof(observed)!="character"){
    stop("ERROR: please enter the column name of the observed ancestry")
  }
  if( !( observed %in% names(data) ) ){
    stop("ERROR: please make sure that the observed ancestry is included in the column names of the data.")
  }
  if( !all(reference %in% names(data) ) ){
    stop("ERROR: Please make sure that all ancestries in reference are also column names in data")
  }
  if(typeof(reference)!="character"){
    stop("ERROR: please enter the column names for the reference ancestries")
  }
  if( !(all(pi.observed<=1) & all(pi.target<=1)) | !(all(pi.observed>=0) & all(pi.target>=0)) ){
    stop("ERROR: pi.observed and pi.target must contain ratios between 0 and 1.")
  }
  k= length(reference)+1
  
  

  #   Normalize pi. We need the pi_reference and pi.target to sum to 1
  pi.target <- pi.target / sum(pi.target)
  pi.observed <- pi.observed / sum(pi.observed)
  
  
  #   We need to sum the reference ancestry allele frequencies multiplied by pi.target. 
  #   We will call this sum "starred"
  #   We also need to sum up the reference ancestry allele frequencies multiplied by pi_reference. 
  #   We will call this sum "hatted"
  
  #   Initialize "hatted" and "starred"
  hatted  <- vector(mode = "double", length = dim(data)[1])
  starred <- vector(mode = "double", length = dim(data)[1])
  
  #   Sum the K-1 reference ancestries multiplied by pi_hat.
  #   Also the sum the k-1 reference ancestries multiplied by pi_star.

  for ( j in 1:(k-1) ){
    hatted  <- hatted  + (pi.observed[j] * data[,reference[j]])
    starred <- starred + (pi.target[j] * data[,reference[j]])
  }
  
  data$adjustedAF <- (
    
    ( pi.target[k] / pi.observed[k] ) *
      
      ( data[,observed] - hatted )
    
    + starred )
  
 ## data.out <- cbind(data[,c("CHR", "RSID", "POS", "REF", "ALT")], data[,c(reference, observed, "adjustedAF")] )
data.out <- data
 
 observed.data=paste("observed data to update AF: '", observed, "'", sep="")
 
 pi_table<-data.frame(ref.group=c(reference, "NONE"), pi.observed=pi.observed, pi.target=pi.target)
 
 Nsnps=nrow(data)
 
 tmp.out<-list("pi"=pi_table, "observed.data"=observed.data, "Nsnps"=Nsnps, "adjusted.AF"=data.out)
 
  print(c(tmp.out[1:3], "use $adjusted.AF to see adjusted AF data"))
 return(tmp.out)

}
