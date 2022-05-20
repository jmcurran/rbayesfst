#' @title Get the current values of the parameters for the priors
#'
#' @param x an object of class bayesFst
#'
#' @return a list containing the following values \tabular{ll}{
#' \code{alpha} \tab a vector containing the mean and std. dev. for the alphas\cr
#' \code{beta} \tab a vector containing the mean and std. dev. for the betas\cr
#' \code{gamma} \tab a vector containing the mean and std. dev. for the gammas\cr
#' \code{uSigma} \tab the scale parameter for the alphas\cr
#' \code{pSigma} \tab the scale parameters for the p's (the allele frequencies)\cr
#' \code{cor} \tab the correlation parameter (fixed) between adjacent loci\cr
#' }
#' @seealso \code{setPriors}
#' @export
getPriors = function(x){
  return(x$b$getPriorParams())
}