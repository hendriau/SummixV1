#' Summix: estimating mixture proportions of reference group
#'
#' Summix: estimating mixture proportions of reference groups from large (N SNPs>10,000) genetic AF data
#'
#' @param data a dataframe of the observed and reference allele frequencies for N genetic variants. See data formatting document at \href{https://github.com/hendriau/Summix}{https://github.com/hendriau/Summix} for more information.
#' @param reference a character vector of the column names for the reference ancestries.
#' @param observed a character value that is the column name for the observed group.
#' @param pi.start length K numeric vector of the starting guess for the ancestry proportions. If not specified, this defaults to 1/K where K is the number of reference ancestry groups.
#' @return data frame with the following columns
#' @return objective: least square value at solution 
#' @return iterations: number of iterations for SLSQP algorithm
#' @return time: time in seconds of SLSQP algorithm
#' @return filtered: number of SNPs not used in estimation due to missing values 
#' @return K columns of mixture proportions of reference ancestry groups input into the function
#'
#' @author Gregory Matesi, \email{gregory.matesi@ucdenver.edu}
#' @author Audrey Hendricks, \email{audrey.hendricks@ucdenver.edu}
#' @reference https://github.com/hendriau/Summix
#' @keywords genetics, mixture distribution, admixture, population stratification
#' 
#' @seealso \code{\link{adjAF}} for adjusting allele frequencies and \url{https://github.com/hendriau/Summix} for further documentation. \code{\link[nloptr]{slsqp}} function in the nloptr package for further details on Sequential Quadratic Programming \url{https://www.rdocumentation.org/packages/nloptr/versions/1.2.2.2/topics/slsqp}

#' 
#' @examples
#' # load the data
#' data("ancestryData")
#' 
#' # Estimate 5 reference ancestry proportion values for the gnomAD African/African American group
#' # using a starting guess of .2 for each ancestry proportion.
#' summix( data = ancestryData, 
#'     reference=c("ref_AF_afr_1000G", 
#'         "ref_AF_eur_1000G", 
#'         "ref_AF_sas_1000G", 
#'         "ref_AF_iam_1000G", 
#'         "ref_AF_eas_1000G"), 
#'     observed="gnomad_AF_afr", 
#'     pi.start = c(.2, .2, .2, .2, .2) )
#'
#' @export
#' @importFrom nloptr slsqp


summix = function(data, reference, observed, pi.start=c()){
  
  if(!is(object = data, class2 = "data.frame")){
    stop("ERROR: data must be a data.frame as described in the vignette")
  }

  if(typeof(observed)!="character"){
    stop("ERROR: 'observed' must be a character string for the column name of the observed ancestry in data")
  }
  if(!(observed %in% names(data))){
    stop("ERROR: 'observed' must be the column name of the observed ancestry in data")
  }
  if(typeof(reference)!="character"){
    stop("ERROR: 'reference' must be a vector of column names from data to be used in the reference")
  }
  if(all(reference %in% names(data))==FALSE){
    stop("ERROR: 'reference' must be a vector of column names from data to be used in the reference")
  }
  
  # Filter NA allele frequencies out of the observed column
  filteredNA <- length(which(is.na(data[,observed]==TRUE)))
  # The math in the ancestr function only uses the 
  # observed allele frequency vector and the reference panel
  observed.b  <- as.data.frame( data[which(is.na(data[,observed])==FALSE),observed] )
  refmatrix <- as.data.frame( data[which(is.na(data[,observed])==FALSE),reference] )
  
  # 
  if(length(pi.start) != 0){
    #########################
    # Check if pi.start is numeric
    if (is.numeric(pi.start)==FALSE){
      stop("ERROR: Please make sure pi.start is a positive numeric vector of length reference that sums to one")
    }
    # Check if length(pi.start)==length(reference)
    if (length(pi.start)!= length(reference)){
      stop("ERROR: Please make sure pi.start is a positive numeric vector of length reference that sums to one")
    }
    # Check that pi.start is positive
    if (all(pi.start>0) == FALSE){
      stop("ERROR: Please make sure pi.start is a positive numeric vector of length reference that sums to one")
    }
    # Check if sum(pi.start) = 1
    if (sum(pi.start)!=1){
      stop("ERROR: Please make sure pi.start is a positive numeric vector of length reference that sums to one")
    }
    #########################
    ###############################
    # Set the starting guess to x_0
    ###############################
    starting = pi.start
  } else{
    starting = rep( 1/ncol(refmatrix), ncol(refmatrix) )
  }
  # Here we are defining the objective function. This function is evaluated at a
  # k-dimensional set  x . Each of our K reference allele frequencies are multiplied by 
  # our current best guess for the ancestry proportion.
  # We then subtract the allele frequency values from the observed homogeneous population.
  # And finally this sum is squared to achieve a least squares form.
  
  #########################
  fn.ancmix = function(x){
    minfunc = 0
    for (i in 1:ncol(refmatrix)){
      minfunc = minfunc + x[i]*refmatrix[,i]
    }
    minfunc = minfunc - observed.b
    minfunc = sum((minfunc)**2)
    return(minfunc)
  }
  ########################
  
  # Here we are defining the gradient of the objective function.
  gr.ancmix <- function(x){
    gradvec = matrix(0,ncol(refmatrix),1)
    gradfunc = 0
    for (i in 1:ncol(refmatrix)){
      gradfunc = gradfunc + x[i]*refmatrix[,i]
    }
    gradfunc = gradfunc - observed.b
    
    for (i in 1:ncol(refmatrix)){
      gradvec[i] = sum(2 * refmatrix[,i] * gradfunc)
    }
    return(gradvec)
  }
  
  # H equality
  # This function returns the equality constraints for the nloptr slsqp algorithm
  # We sum up the K current proportion estimate values and subtract 1. If the estimated proportion values sum to 1 than this value should equal zero.
  heq.ancmix = function(x){
    equality = 0
    for (i in 1:ncol(refmatrix)){
      equality = equality + x[i]
    }
    return(equality - 1)
  }
  
  # H inequality
  # This function returns a K vector of 
  hin.ancmix <- function(x){
    h = numeric(ncol(refmatrix))
    for (i in 1:ncol(refmatrix)){
      h[i] = x[i]
    }
    return(h)
  }
  
  # We use the start_time function base function 
  #to record the run time for our convex optimization algorithm.
  # The output for the nloptr slsqp function is stored in the variable S.
  # This function requires 5 inputs:
  #   1. fn.ancmix is the objective function.
  #   2. gr.ancmix is the gradient of the objective function.
  #   3. hin.ancmix defining the inequality constraints
  #   4. heq.ancmix defining the equality constraints
  start_time = Sys.time()
  S = suppressMessages( slsqp(starting,
            fn = fn.ancmix,
            gr = gr.ancmix,
            hin = hin.ancmix,
            heq = heq.ancmix)
      )
  end_time = Sys.time()
  ttime = end_time - start_time
  
  # par is the K optimal ancestry propotion estimates
  # value is the minimization function evaluated at par
  # iter is the number of iterations that the algorithm took to reach the optimal solution of par
  # finally ttime is the run time for the algorithm
  
  d <- data.frame(matrix(ncol = length(reference)+4, nrow = 1)) 
  colnames(d) <- c("objective", "iterations", "time",
                        "filtered", colnames(refmatrix))
  
  d[1] <- S$value
  d[2] <- S$iter
  d[3] <- ttime
  d[4] <- filteredNA
  for(i in 1:length(reference)){
    d[4+i] <- S$par[i]
  }
  return(d)
}
