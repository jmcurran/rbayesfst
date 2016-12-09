readData = function(fileName){
  f1 = file(fileName, "r")

  if(!isOpen(f1)){
    stop(paste0("Cannot read ", fileName))
  }

  Lines = readLines(f1)
  close(f1)

  gammaSwitch = as.logical(Lines[1])
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

  return(list(gammaSwitch = gammaSwitch, nPops = nPops, nLoci = nLoci, Loci = Loci, Pops = Pops,
              numAlleles = numAlleles, dbCounts = dbCounts, locusPopSums = locusPopSums))
}
