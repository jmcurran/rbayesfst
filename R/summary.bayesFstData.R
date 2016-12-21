#' @title Summary information for \code{bayesFstData} objects
#' 
#' @description Provides a simple summary of the input data used with this package.
#' 
#' @param x An object of class \code{bayesFstData}.
#' @param \ldots Not used.
#' 
#' @export
#' @seealso readData
summary.bayesFstData = function(x, ...){
  cat(sprintf("This dataset contains %d populations and %d loci\n", x$nPops, x$nLoci))
  cat("The populations in this dataset are:\n")
  print(x$Pops)
  cat("\n\nThe loci in this dataset are:\n")
  print(x$Loci)
  cat(sprintf("The minimum number of alleles per locus is %d\n", min(x$numAlleles)))
  cat(sprintf("The maximum number of alleles per locus is %d\n", max(x$numAlleles)))
}