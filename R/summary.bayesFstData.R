#' An S3 summary method for bayesFstData objects
#' 
#' Provides a simple summary of the input data used with this package.
#' 
#' @param x An object of class \code{bayesFstData}.
#' @param \ldots Not used.
#' 
#' @export
summary.bayesFstData = function(x, ...){
  cat(sprintf("This dataset contains %d populations and %d loci\n", x$nPops, x$nLoci))
  cat("The populations in this dataset are:\n")
  print(x$Pops)
  cat("\n\nThe locus in this dataset are:\n")
  print(x$Loci)
  cat(sprintf("The minimum number of alleles per locus is %d\n", min(l$numAlleles)))
  cat(sprintf("The maximum number of alleles per locus is %d\n", max(l$numAlleles)))
}