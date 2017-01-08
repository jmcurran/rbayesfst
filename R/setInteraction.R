#' @title Change the model to use (or not use) the interaction term
#'
#' @param x an object of class \code{bayesFst}
#' @param bInteraction if \code{TRUE} then the interaction (gamma) terms are used in the model, otherwise they are not.
#' 
#' @note This function changes a value in a C++ object. As a consequence it DOES NOT need to be reassigned.
#' @export
#'
#' @examples
#' bd = readData(system.file("extdata", "data_BB04.json", package = "rbayesfst"))
#' bf = init(bd)
#' print(bf)
#' ## turn the interaction term on
#' setInteraction(bf)
#' ## note how the model has changed even though there is no reassignment
#' print(bf)
setInteraction = function(x, bInteraction = TRUE){
  x$b$interaction = bInteraction
}