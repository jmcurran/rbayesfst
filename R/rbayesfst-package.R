#' rbayesfst
#' 
#' An R interface to the BayesFst programme. This is a Metropolis algorithm for estimating F_ST, the correlation of two homologous genes sampled within the same subpopulation, relative to the total population. Among other uses, F_ST is important in calculating forensic match probabilities, as matches may be more common between individuals in the same subpopulation (due to recent shared ancestry of genes in the DNA profile), see Balding & Nichols (1994, Forensic Science International 64, pp125-140). F_ST is also interesting as an indicator of possible selection, see Lewontin & Krakauer (1973, Geneticsc74:175-195) and citations thereof, especially Beaumont & Balding (2004, Molecular Ecology, 13:969-980).
#'
#' @docType package
#' @name rbayesfst
#' @title Bayesian Estimation of Fst
#' 
#' @author James Curran, David Balding, Mark Beaumont
#' @useDynLib rbayesfst, .registration = true
NULL
