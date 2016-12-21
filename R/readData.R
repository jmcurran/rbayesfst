#' readData
#' 
#' @description A function to read the data needed for the algorithm.
#' 
#' @details A function that can read the data needed for the algorithm. The original program requires data to be in the following format 
#' Input format
#' 
#' The input file should be a space or tab-delimited ASCII (i.e. text)
#' file. First, there are three header lines, containing:
#' 
#' \enumerate{
#' \item interaction term (gamma) switch: 0 - for no interaction term 1 - to include interaction term
#' \item the number of subpopulations (npop)
#' \item the number of loci (nloc)
#' }
#' 
#' Then for locus i (i=1,..,nloc) there are npop+1 input lines, the first
#' containing the number of alleles at that locus (numall[i]) and the
#' remaining lines each contain the allele counts.
#' 
#' For example, for 3 subpops and 2 loci, with 2 and 4 alleles, where we
#' do not want to include the interaction gamma in the analysis, the input
#' should look like:
#'   
#' \code{0} \cr
#' \code{3} \cr
#' \code{2} \cr
#' 
#' 
#' \code{2} \cr
#' \code{45 5} \cr
#' \code{23 27} \cr
#' \code{10 40} \cr
#'  \cr
#' \code{4} \cr
#' \code{10 10 10 20} \cr
#' \code{5 5 25 15} \cr
#' \code{17 19 0 14} \cr \cr
#  
#' Blank lines are used to delimit information for the different loci. White space (space or tab) is used to delimit the numbers on each line.
#' The locus and population names are extracted from comments on line which start with a \code{\%}, e.g. using the same example as above
#' 
#' \code{0} \cr
#' \code{3} \cr
#' \code{2} \cr
#' 
#' 
#' \code{2} \cr
#' \code{45 5 \% Locus1 Pop1}\cr
#' \code{23 27 \% Locus1 Pop2}\cr
#' \code{10 40 \% Locus1 Pop3}\cr
#'  \cr
#' \code{4} \cr
#' \code{10 10 10 20 \% Locus2 Pop1}\cr
#' \code{5 5 25 15 \% Locus2 Pop2}\cr
#' \code{17 19 0 14 \% Locus2 Pop3}\cr \cr
#'
#' The function also has the ability to read in JSON formatted data. As an absolute minimum it must contain the the counts of each allele at each locus for each population stored in a \code{list} called
#' \code{dbCounts}. Each element of the list should be a \code{matrix} with the same number of rows, one for each population. The number of columns in each row is equal to the number of alleles for the 
#' particular locus, and the entries of each row give the count of allele \eqn{i}{i} in population \eqn{j}{j} at locus \eqn{l}{l}.
#' 
#' @param fileName A valid file path or URL with the file extension \code{inp} or \code{json}.
#'
#' @return A dataset of class \code{bayesFst}. This is essentially a list with the following members: \tabular{ll}{
#' \code{dbCounts} \tab A list of matrices providing allele counts at each locus for each population. The items in the list are loci, the rows are populations, and the columns alleles. \cr
#' \code{nLoci} \tab the number of loci in the dataset. \cr
#' \code{Loci} \tab the locus names if available. If the locus names are not available they will be labelled Locus 1, Locus 2, \ldots \cr
#' \code{nPops} \tab the number of populations in the dataset. \cr
#' \code{Pops} \tab the population names if available. If the population names are not available they will be labelled Pop 1, Pop 2, \ldots \cr
#' \code{numAlleles} \tab the number of possible alleles at each locus \cr
#' \code{locusPopSums} \tab a \eqn{n_{loc}\times n_{pop}}{nloc x npop} matrix containing the number of alleles observed at the \eqn{l^{\mathrm{th}}}{lth} locus for the \eqn{j^{\mathrm{th}}}{jth} population. \cr
#' \code{gammaSwitch} \tab either \code{TRUE} or \code{FALSE} depending on whether locus and popuation effects interact or not.\cr
#' }
#' @export
#' 
#' @examples 
#' \dontrun{
#' bd = readData('http://www.reading.ac.uk/Statistics/genetics/software/bayesfst/data_BB04.inp')}
#' 
#' @seealso summary.bayesFstData
readData = function(fileName){
  if(grepl("^.*inp$", fileName)){
    cat("Reading raw data in .inp format\n")
    
    f1 = file(fileName, "r")
  
    if(!isOpen(f1)){
      stop(paste0("Cannot read ", fileName))
    }
  
    Lines = readLines(f1)
    close(f1)
  
    gammaSwitch = Lines[1] == 1
    nPops = as.numeric(Lines[2])
    nLoci = as.numeric(Lines[3])
    Lines = Lines[-c(1:3)]
  
    ## blank line means new locus
    idx = grep("^$", Lines)
    if(length(idx) != nLoci){
      stop("Number of loci does not match")
    }
  
    numAlleles = rep(0, nLoci)
    Loci = rep("", nLoci)
    Pops = rep("", nPops)
    dbCounts = vector(length = nLoci, mode = "list")
  
    for(loc in 1:nLoci){
      start = idx[loc] + 1
      end = if(loc < nLoci){ idx[loc + 1] - 1} else {length(Lines)}
      locusLines = Lines[start:end]
  
      numAlleles[loc] = as.numeric(locusLines[1])
      locusLines = locusLines[-1]
      dbCounts[[loc]] = matrix(0, nrow = nPops, ncol = numAlleles[loc])
      
      for(pop in 1:nPops){
        line = locusLines[pop]
        splitLine = unlist(strsplit(line, "%"))
        counts = as.numeric(unlist(strsplit(gsub('[[:space:]]+', ',',  splitLine[1]), ',')))
        dbCounts[[loc]][pop, ] = counts
  
        locusPopDetails = unlist(strsplit(gsub("^[[:space:]]+([a-zA-Z0-9\\-]+)[[:space:]]+([A-Z_]+)[[:space:]]*$", "\\1,\\2", splitLine[2]), ","))
  
        if(pop == 1){
          Loci[loc] = locusPopDetails[1]
        }
  
        if(loc == 1){
          Pops[pop] = locusPopDetails[2]
        }
      }
    }
  

    locusPopSums = do.call(rbind, lapply(dbCounts, function(x)rowSums(x)))
    l = list(gammaSwitch = gammaSwitch, nPops = nPops, nLoci = nLoci, Loci = Loci, Pops = Pops,
         numAlleles = numAlleles, dbCounts = dbCounts, locusPopSums = locusPopSums)
    class(l) = "bayesFstData"

    return(l)
  }else if(grepl("^.*json$", fileName, ignore.case = TRUE)){
    cat("Reading file in JSON format\n")
    l = fromJSON(readLines(fileName))
    
    nms = names(l)
    
    if(!("dbCounts" %in% nms))
      stop("This data set doesn't contain any allele counts. This is the minimal set of information needed for this programme.")
    
    if(!("Loci" %in% nms)){
      l$nLoci = length(l$dbCounts)
      
      if(is.null(names(l$dbCounts))){
        l$Loci = paste("Locus", 1:(l$nLoci))
      }else{
        l$Loci = names(l$dbCounts)
      }
    }
    
    if(!("Pops" %in% nms)){
      l$nPops = nrow(l$dbCounts[[1]])
      
      if(is.null(rownames(l$dbCounts[[1]]))){
        l$Pops = paste("Pop", 1:(l$nPops))
      }else{
        l$Pops = rownames(l$dbCounts[[1]])
      }
    }
    
    if(!("numAlleles" %in% nms)){
      l$numAlleles = sapply(l$dbCounts, ncol)
    }
    
    if(!("locusPopSums" %in% nms)){
      l$locusPopSums = do.call(rbind, lapply(l$dbCounts, function(x)rowSums(x)))
    }
    
    class(l) = "bayesFstData"
    
    return(l)
  }else{
    stop("I don't know how to read this file format.\n")
  }
}
