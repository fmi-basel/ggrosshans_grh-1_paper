#
# Helper functions for assigning regions to regions
#


#' Wrapper designed to assign one set of query genomic regions to another set of
#' subject genomic regions.
#' Assignment is based on region overlaps, or preceding/following the subjects.
#' Then calculates overlap length and distance between assigned regions;
#' supports multiple ways to calculate said distance.
#' 
#' Typical use case would be using ChIP-seq peaks as queries and gene promoter regions/
#' bodies as subjects, using a summit to TSS distance.
#' 
#' @param queries Named GRanges object.
#' @param subjects Named GRanges object.
#' @param assign_type Type of assignment (overlap, precede, or follow).
#' @param query_IDs Variable name in which to store names of query GRanges entries.
#' @param subject_IDs Variable name in which to store names of subject GRanges entries.
#' 
#' @param query_dist_type,subject_dist_type Determines how distances between ranges
#'  in queries and subjects are calculated. 
#'  "range" (default) simply uses he full range.
#'  "tss" takes the earliest nucleotide in the range (eg where the range ends if strand is -).
#'  "mid" takes the midpoint of the range.
#'  "summit" looks for a "summit" column in the GRanges object, and considers that
#'  to be an offset from the range start, and takes the resulting nucleotide.
#'  
#' @param verbose Boolean, should some progress be printed to screen?
#' 
#' @return An annotated Hits object.
#' 
#' @seealso nearest-methods of package GenomicRanges.
assign_ranges <- function(
  
  queries,
  subjects, 
  assign_type = c("overlap", "precede", "follow"),
  
  query_IDs="peakName", subject_IDs="geneID",
  
  query_dist_type   = c("range", "tss", "start", "center", "end", "mid", "summit"),
  subject_dist_type = c("range", "tss", "start", "center", "end", "mid", "summit"),
  
  verbose = F, 
  
  select = "all",  # passed to findOverlaps() or precede()
  ...              # passed to findOverlaps() or precede()
  
  ) {
  
  # checks
  assign_type <- match.arg(assign_type) 
  query_dist_type <- match.arg(query_dist_type) 
  subject_dist_type <- match.arg(subject_dist_type) 
  # queries & subjects MUST be *named* GRanges objects
  stopifnot(inherits(queries, "GRanges"))
  stopifnot(inherits(subjects, "GRanges"))
  stopifnot(all( c( !is.null(names(queries)),
                    !is.null(names(subjects))  )))
  
  if (verbose) {
    cat("\nSUBJECTS:\n")
    print(subjects)
    cat("\nQUERIES:\n")
    print(queries)
  }
  
  
  # find overlaps
  assign_hits <- switch(
    assign_type,
    "overlap" = findOverlaps( 
      queries, subjects, select=select, ...),
    "precede" =  { 
      stopifnot(select=="all") 
      precede(     
        queries, subjects, select=select, ...)
      },
    "follow"  =  follow(     
      queries, subjects, select=select, ...)
  )
  # assign_hits is an object inheriting from the Hits class
  
  if (length(assign_hits)<=0) { return(assign_hits) }
  
  
  # transform indices to IDs
  mcols(assign_hits)[[query_IDs]]   <- 
    names(queries)[  queryHits(assign_hits)   ]
  mcols(assign_hits)[[subject_IDs]] <- 
    names(subjects)[ subjectHits(assign_hits) ]
  
  
  # record assignment type
  mcols(assign_hits)[,"assignmentProcess"] <- paste(
    query_dist_type, assign_type, subject_dist_type, sep="_") 
  
  
  # calculate overlap lengths
  mcols(assign_hits)[["overlap"]] <- switch(
    assign_type,
    "overlap" = width( pintersect(
      queries[from(assign_hits),], subjects[to(assign_hits),]) ),
    "precede" =  0,
    "follow" =  0
  )
  
  
  # establish ranges to use for distance calculations
  get_dist_granges <- function(mygr, opt) {
    if (opt %in% c("start", "center", "end")) {
      resize(mygr, width=1, fix=opt)
    } else {
      switch(
        opt,
        "range" = mygr,
        "tss"   = {
          warning("Assuming that range strand-sensitive start is TSS.\n")
          resize(mygr, width=1, fix="start")
        },
        "mid"   = resize(mygr, width=1, fix="center"),
        "summit" = {
          stopifnot( "summit" %in% colnames(mcols(mygr)))
          start(mygr) <-
            start(mygr) + mcols(mygr)[,"summit"] # NB: assumes to be relative
          resize(mygr, width=1, fix="start")
        }
      )
    }
  }
  subject_dist_gr <- get_dist_granges(subjects, subject_dist_type)
  query_dist_gr   <- get_dist_granges(queries,  query_dist_type)

  if (verbose) {
    cat("\nSUBJECT DISTS:\n")
    print(subject_dist_gr)
    cat("\nQUERY DISTS:\n")
    print(query_dist_gr)
  }


  # calculate distances
  # direction based on distance GR midpoints, negative if query position is lower than subject
  mcols(assign_hits)[,"distance"] <- distance(
    query_dist_gr[from(assign_hits),],
    subject_dist_gr[to(assign_hits),])
  bNeg <- mid(query_dist_gr[from(assign_hits),]) < mid(subject_dist_gr[to(assign_hits),])
  mcols(assign_hits)[bNeg,"distance"] <- -1 * mcols(assign_hits)[bNeg,"distance"]

  # return
  if (verbose) {
    cat("\nRESULTS:\n")
    print(assign_hits)
  }
  return(invisible(assign_hits))
  
}


if (FALSE) {
  ## overlap test ##
  gr_s <- GRanges(
    Rle(c("chr1"), c(1)),
    IRanges(20, width=50),
    strand=Rle(strand(c("-")), c(1)) )
  names(gr_s) <- paste0("subj",1:length(gr_s))
  
  gr_q <- GRanges(
    Rle(c("chr1"), c(1)),
    IRanges(15, width=10),
    strand="*")
  gr_q$summit <- 2
  names(gr_q) <- paste0("quer",1:length(gr_q))
  
  assign_ranges(
    gr_q,
    gr_s,
    overlap_type = "overlap",
    query_dist_type = "summit",
    subject_dist_type = "tss"
  )
}

assign_ranges_QC <- function(
    assign_hits, 
    queries, subjects, 
    main_blurb="", dist_blurb="") {
  
  stopifnot(inherits(assign_hits, "Hits"))
  
  mcols(assign_hits)$strand <- strand(subjects[mcols(assign_hits)$geneID])
  
  h <- function(strd) {
    ahs <- subset(assign_hits, strand==strd)
    hist(mcols(ahs)[,"distance"], breaks=32, 
         main=paste0(main_blurb, " (",strd," strand)"),
         xlab=paste0(dist_blurb, " distance"),
         sub=paste(
           length(unique(from(ahs))), "/", length(queries),
           "peaks assigned to",
           length(unique(to(ahs))), "/", length(subjects),
           "genes" ) )
    abline(v = 0, col="red", lty="dashed", lwd=3)
  }
  h( "+" )
  h( "-" )
  
  if (any(na.omit(mcols(assign_hits)[["overlap"]])>0)) {
    hist(mcols(assign_hits)[["overlap"]], breaks=64, 
         main=paste0(main_blurb),
         xlab="overlap")
    abline(v = 0, col="grey10", lty="dashed", lwd=3)
  }
  
}


if (FALSE) {
  ## preceding test ##
  gr_s <- GRanges(
    Rle(c("chr1"), c(1)),
    IRanges(30, width=50),
    strand=Rle(strand(c("+")), c(1)) )
  names(gr_s) <- paste0("subj",1:length(gr_s))
  
  gr_q <- GRanges(
    Rle(c("chr1"), c(1)),
    IRanges(15, width=10),
    strand="*")
  gr_q$summit <- 2
  names(gr_q) <- paste0("quer",1:length(gr_q))
  
  assign_ranges(
    gr_q,
    gr_s,
    overlap_type = "precede",
    query_dist_type = "summit",
    subject_dist_type = "tss"
  )
}

collate_assignments <- function(
    ar_res_list, 
    ID_cols=c("peakName","geneID") ) {
  
  fmt_hits <- function(h, type) {
    as.data.frame(cbind(
      type=type,
      mcols(h)[,c(ID_cols,"overlap","distance")]
    ))
  }
  
  for (n in names(ar_res_list)) {
    ar_res_list[[n]] <- fmt_hits(ar_res_list[[n]], n)
  }
  
  res <- as.data.frame(do.call(rbind, ar_res_list))
  rownames(res) <- NULL
  return(res)
}






# EOF
