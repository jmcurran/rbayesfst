#' @title Summary information for \code{bayesFstData} objects
#' 
#' @description Provides a simple summary of the input data used with this package.
#' 
#' @param object An object of class \code{bayesFstData}.
#' @param \ldots Not used.
#' 
#' @export
#' @seealso readData
summary.bayesFstData = function(object, ...){
  cat(sprintf("This data is loaded from the file: \n%s\n\n", object$name))
  cat(sprintf("This dataset contains %d populations and %d loci\n", object$nPops, object$nLoci))
  cat("The populations in this dataset are:\n")
  print(object$Pops)
  cat("\n\nThe loci in this dataset are:\n")
  print(object$Loci)
  cat(sprintf("The minimum number of alleles per locus is %d\n", min(object$numAlleles)))
  cat(sprintf("The maximum number of alleles per locus is %d\n", max(object$numAlleles)))
}