#' sample generic
#' 
#' @param x an object.
#' @param \dots Any additional arguments to be passed to \code{sample}.
#' @export
sample = function(x, ...){
  UseMethod("sample")
}
#' Default call to base R \code{\link[base]{sample}} function
#' 
#' @param x either a vector of one or more elements from which to choose,
#' or a positive integer. 
#' @param size a non-negative integer giving the number of items to choose.
#' @param replace should sampling be with replacement?
#' @param prob a vector of probability weights for obtaining the elements of 
#' the vector being sampled.
#' @param \ldots any other arguments. Not used.
#' @seealso \code{\link[base]{sample}}
#' @export
sample.default = function(x, size, replace = FALSE, prob = NULL, ...){
  base::sample(x, size, replace, prob)
}
#' Sample from the posterior distribution of Fst
#' @param x an object class \code{bayesFst}
#' @param size the desired numbers of samples to be taken from the posterior distribution of the parameters
#' @param seed a positive (32 bit) integer to be used as a random number seed in a Mersenne-Twister random number generator. This package uses the C++ `random` library.
#' @param \ldots not used
#' @note If you want to change some of the parameters for the MCMC algorithm, other than the number of outputs, then use the specific methods to do this. They cannot be changed here.
#' @seealso \code{\link{setInteraction}}, \code{\link{setPriors}}, \code{\link{setRunParams}}
#' @importFrom coda mcmc
#' @importFrom stats runif
#' @return a list 
#' @export
sample.bayesFst = function(x, size = NULL, seed = floor(runif(1, 0, 2^32) + 1), ...){
  if(!is.null(size)){
    setRunParams(x, numOut = size)
  }
  results = x$b$run(seed)
  
  thin = x$b$getThin()
  
  colnames(results$alpha) = x$data$Loci
  results$alpha = mcmc(results$alpha, thin = thin)
  
  colnames(results$beta) = x$data$Pops
  results$beta = mcmc(results$beta, thin = thin)
  
  if(x$b$interaction){
    colnames(results$gamma) = paste((rep(x$data$Loci, rep(length(x$data$Pops), 
                                                         length(x$data$Loci)))), x$data$Pops, sep = "x")
    results$gamma = mcmc(results$gamma, thin = thin)
  }
  return(results)
}