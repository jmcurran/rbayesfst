#' @title Set the values of the parameters for the priors
#'
#' @param x an object of type bayesFst
#' @param alphaMu (expected) mean of Normal distribution for prior of the locus effect alpha[i]
#' @param alphaSigma (expected) mean of Normal distribution for prior of locus effect alpha[i]
#' @param betaMu mean of Normal distribution for prior of the population effect beta[j]
#' @param betaSigma sd of Normal distribution for prior of the population effect  beta[j]
#' @param gammaMu  mean of Normal distribution for prior of the locus x population interaction effect gamma[i,j]
#' @param gammaSigam  sd of Normal distribution for prior the locus x population interaction effect of gamma[i,j]
#' @param uSigma scale param for Normal updates
#' @param pSigma scale param for Dirichlet updates of the allele frequencies p[i,j]
#' @param cor The correlation parameter (fixed) between adjacent loci
#' 
#' @note This function directly affects a C++ object and therefore DOES NOT return anything.
#' @seealso getPriors
#'
#' @export
#'
#' @examples
#' bd = readData(system.file("extdata", "data_BB04.json", package = "rbayesfst"))
#' bf = init(bd)
#' print(bf)
#' ## change the alpha mean only
#' setPriors(bf, alphaMu = 1)
#' ## note how the mean has changed even though there is no reassignment
#' print(bf)
setPriors = function(x, alphaMu = 0.0, alphaSigma = 1.0,  betaMu = -2.0, betaSigma = 1.8, 
                     gammaMu = 0.0, gammaSigma = 0.5, uSigma = 0.5, pSigma = 1000, cor = 0){
  x$b$setPriorParams(c(alphaMu, alphaSigma), c(betaMu, betaSigma), c(gammaMu, gammaSigma),
                         uSigma, pSigma, cor)
}