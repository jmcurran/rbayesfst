#' sample generic
#' 
#' @param x an object.
#' @param \dots Any additional arguments to be passed to \code{sample}.
#' @export
sample = function(x, ...){
  UseMethod("sample")
}
#' Default call to base R \code{\link[base]{sample}} function
#' @seealso \code{\link[base]{sample}}
#' @export
sample.default = function(x, size, replace = FALSE, prob = NULL){
  base::sample(x, size, replace, prob)
}
#' Sample from the posterior distribution of Fst
#' @param x an object class \code{bayesFst}
#' @param \ldots not used
#' @note If you want to change some of the parameters for the MCMC algorithm, other than the number of outputs, then use the specific methods to do this. They cannot be changed here.
#' @seealso \code{\link{setInteraction}}, \code{\link{setPriors}}, \code{\link{setRunParams}}
#'  @return a list 
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