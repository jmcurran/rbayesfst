#' Intialise a new Bayes Fst Object
#'
#' @param data a dataset of class \code{bayesFstData} loaded by \code{readData}.
#'
#' @return an object of class \code{bayesFst}
#' 
#' @examples 
#' bd = readData(system.file("extdata", "data_BB04.json", package = "rbayesfst"))
#' bf = init(bd)
#' @export
init = function(data){
  b = new(BayesFst)
  b$setData(data)
  
  myObj = list(b = b)
  
  class(myObj) = "bayesFst"
  
  return(myObj)
}