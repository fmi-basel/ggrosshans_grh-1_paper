#

# drop-in replacement for monaLisa::homerToPFMatrixList
# that can't handle invalid matrices sometimes returned by HOMER,
# and that can't pass other things to TFBSTools::PFMatrix,
# eg BG frequencies.
homerToPFMatrixList <- function (
    filename, 
    n=1000, # as smallest frequency from HOMER is 0.001, using n=1000 makes sense
    ... ) {
  
  # set up
  monaLisa:::.assertScalar(x = filename, type = "character")
  stopifnot(file.exists(filename))
  monaLisa:::.assertScalar(x = n, type = "numeric", rngExcl = c(0, Inf))
  read_lines <- readLines(filename)
  
  # ID header lines
  g <- grep(">", read_lines)
  read_lines[g] <- sub("^>", "", read_lines[g])
  fields <- strsplit(read_lines[g], "\t", fixed = TRUE)
  
  # iterate across headers
  L <- lapply(seq_along(g), function(i) {
    
    # extract motif info
    # CF http://homer.ucsd.edu/homer/motif/index.html
    cons <- fields[[i]][1]
    nm <- fields[[i]][2]
    log2cut <- round(as.numeric(fields[[i]][3])/log(2), 2)
    pval <- exp(as.numeric(fields[[i]][4]))
    
    # locate & extract matrix values
    s <- g[i] + 1L
    e <- if (i == length(g)) 
      length(read_lines)
    else g[i + 1L] - 1L
    m <- round(n * do.call(cbind, lapply(strsplit(
      read_lines[s:e], "\t", fixed = TRUE), as.numeric)), 0)
    rownames(m) <- c("A", "C", "G", "T")
    if(ncol(m)==13) 
      print(m)
    
    # assemble PFM
    if (any(is.na(m))) {
      warning("HOMER motif PWM found containing invalid values: ",nm)
      return(NULL)
    }
    TFBSTools::PFMatrix(
      ID = nm, name = nm, profileMatrix = m, 
      tags = list(
        log2cut = log2cut, 
        pval    = pval,
        comment = "imported from HOMER motif file"),
      ... )
  })
  L <- L[!unlist(lapply(L,is.null))]
  
  # return
  do.call(PFMatrixList, c(L, list(use.names = TRUE)))
}


# EOF