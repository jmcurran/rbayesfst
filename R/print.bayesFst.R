#' @title Print method for \code{bayesFst} objects
#' 
#' @description Provides a simple summary of the input data used with this package.
#' 
#' @param x An object of class \code{bayesFstData}.
#' @param \ldots Not used.
#' 
#' @export
#' @seealso readData
print.bayesFst = function(x, ...){
  cat("This object provides the input parameters and data for the Bayesian estimation\n")
  cat("of Fst. The model is\n\n")
  cat("\tlogit(Fst[i,j]) = alpha[i] + beta[j]")
  
  if(x$b$interaction){
    cat(" + gamma[i,j]\n\n")
  }else{
    cat("\n\n")
  }
  
  cat("where:\n")
  cat("\talpha[i] is the effect for locus i,\n")
  cat("\tbeta[j] is the effect for population j")
  
  
  if(x$b$interaction){
    cat(",\n\tgamma[i,j] is the effect for locus i in population j\n")
  }else{
    cat(".\n")
  }
  
  cat("\nThe current priors are:\n\n")
  l = x$b$getPriorParams()
  cat(sprintf("  alpha[i] ~ N(%.2f, %.2f^2),\n", l$alpha[1], l$alpha[2]))
  cat(sprintf("   beta[i] ~ N(%.2f, %.2f^2)", l$beta[1], l$beta[2]))
  
  if(x$b$interaction){
    cat(sprintf(",\ngamma[i,j] ~ N(%.2f, %.2f^2),\n\n", l$gamma[1], l$gamma[2]))
  }else{
    cat(".\n\n")
  }
  
  cat(sprintf("The scale parameter for Normal updates is %.2f\n", l$uSigma))
  cat(sprintf("The scale parameter for Dirichlet updates of p[i,j] is %.2f\n", l$pSigma))
  cat(sprintf("The correlation parameter (fixed) between adjacent loci is %.2f\n\n", l$cor))
  
  x$b$printRunInfo()
  
  cat("\nData has ")
  if(!x$b$isDataLoaded()){
    cat("not ")
  }
  cat("been loaded.\n")
}