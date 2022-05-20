#' @title Set the parameters for the MCMC run
#'
#' @param x an object of class bayesFst
#' @param numOut the desired numbers of samples to be taken from the posterior
#' @param keepPpn \eqn{\pi_k} a thinning proportion used to calculate the interval between successive outputs \eqn{k = \pi_k(n_P+n_L)} where \eqn{n_P} is the number of populations and \eqn{n_L} is the number of loci. This number is in turn used to calculate the number of interations after burn-in \eqn{n_{iter}=k\times n_out}
#' @param discardPpn \eqn{\pi_{discard}} the proportion of iterations to be used as burn-in used to calculate \eqn{d = n_{iter}\pi_d}. The total number of iterations is \eqn{d+n_{iter}}
#' @param acceptPpn \eqn{\pi_a} used to calculate the gap between output of acceptance rates \eqn{a = \pi_a(d+n_iter)}
#' @param print controls the output of something but I forget what at this point
#'
#' @note This function sets values in a C++ object. As a consequence it DOES NOT need to be re-assigned.
#' @export
#'
#' @examples
#' bd = readData(system.file("extdata", "data_BB04.json", package = "rbayesfst"))
#' bf = init(bd)
#' ## set things up for a short test run
#' setRunParameters(bd, numOut = 20)
setRunParams = function(x, 
                        numOut = 2001, 
                        keepPpn = 0.2, 
                        discardPpn = 0.05, 
                        acceptPpn = 0.02, 
                        print = FALSE){
  x$b$setRunParameters(numOut, keepPpn, discardPpn, acceptPpn, print)
}