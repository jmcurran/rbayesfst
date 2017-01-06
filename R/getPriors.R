#' @title Get the current values of the parameters for the priors
#'
#' @param x an object of class bayesFst
#'
#' @return a list containing the following values \tabular{ll}{
#' \code{alpha} \tab a vector containing the mean and std. dev. for the alphas\\
#' \code{beta} \tab a vector containing the mean and std. dev. for the betas\\
#' \code{gamma} \tab a vector containing the mean and std. dev. for the gammas\\
#' \code{uSigma} \tab the scale parameter for the alphas\\
#' \code{pSigma} \tab the scale parameters for the p's (the allele frequencies)\\n
#' \code{cor} \tab the correlation parameter (fixed) between adjacent loci
#' }
#' @seealso \code{setPriors}
#' @export
getPriors = function(x){
  return(x$b$getPriorParams())
}