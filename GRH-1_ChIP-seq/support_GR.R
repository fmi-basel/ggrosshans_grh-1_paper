#
# Support functions for common manipulations on GRanges objects
#

library(stringr)

drop_mcols <- function(gr) {
  mcols(gr) <- NULL
  gr
}


load_MACS2_to_GR <- function(fpn, genome=BSgenome.Celegans.UCSC.ce11) {
  
  peak <- data.table::fread(
    fpn, header = F, data.table = F )
  colnames(peak) <- c(
    "chrom", "start", "end", "name", "score", "strand", 
    "MACS2_foldEnr","log10pval","MACS2_log10qVal","summit")
  if (!all(str_detect(peak$chrom, "chr")))
    peak$chrom <- paste0("chr",peak$chrom)
  peak <- subset(peak, chrom %in% seqnames(genome))
  
  peaks_gr <- GenomicRanges::makeGRangesFromDataFrame(
    peak, keep.extra.columns = T, seqinfo = seqinfo(genome))
  names(peaks_gr) <- peaks_gr$name
  
  return(peaks_gr)
}

getGC <- function( gr, genome=BSgenome.Celegans.UCSC.ce11 ) {
  gr_seq <- getSeq(genome, gr)
  nt_cnts <- alphabetFrequency(gr_seq)[,c("A","C","T","G")]
  gr$pGC <- 100*rowSums(nt_cnts[,c("C","G")]) / rowSums(nt_cnts)
  return(gr) 
}

getMap <- function(gr, genome="ce11") {
  
  # Prepared by Michael using a Swissknife function:
  #   https://fmicompbio.github.io/swissknife/reference/getMappableRegions.html
  # 
  # Essentially, any 50-mer in the genome that does not map to a unique location is
  # called as unmappable. The whole genome has 1bp resolution for this, and the results
  # were processed to generate a GRanges object of contiguous mappable bases.
  # 
  # You can overlap these regions to derive %mappability for any given query region.
  
  mapgr <- FMIRegionDB::getMappableRegions(genome)
  strand(mapgr) <- "*" # was prepared stranded but doesn't really matter
  
  hits <- findOverlaps(gr, mapgr)
  mcols(hits)$overlap <- width( pintersect(
    gr[from(hits),], mapgr[to(hits),]) )
  hits <- as.data.frame(hits)
  tmp <- aggregate(hits$overlap, hits[,"queryHits",drop=F], sum)
  
  gr$pMap <- NA
  mcols(gr)[tmp$queryHits,"pMap"] <- tmp$x
  gr$pMap <- 100 * gr$pMap / width(gr)
  gr$pMap[is.na(gr$pMap)] <- 0
  return(gr)
  
}


prep_nt_table <- function(obj, genome=BSgenome.Celegans.UCSC.ce11, asfactor=F) {
  
  if (inherits(obj, "GRanges")) {
    seqs <- getSeq(genome, obj)
    nt_df <- as.data.frame(do.call(rbind, str_split(seqs, "")))
    colnames(nt_df) <- paste0("pos",1:ncol(nt_df))
    rownames(nt_df) <- names(seqs)
  } else if (inherits(obj, "MsaDNAMultipleAlignment")) {
    seqs <- as.character(obj)
    nt_df <- as.data.frame(do.call(rbind, str_split(seqs, "")))
    colnames(nt_df) <- paste0("pos",1:ncol(nt_df))
    rownames(nt_df) <- names(seqs)
  } else if (inherits(obj, "character")) {
    seqs <- obj
    nt_df <- as.data.frame(do.call(rbind, str_split(seqs, "")))
    colnames(nt_df) <- paste0("pos",1:ncol(nt_df))
    rownames(nt_df) <- names(seqs)
  }
  
  if (asfactor)
    for (i in 1:ncol(nt_df))
      nt_df[,i] <- factor(nt_df[,i], levels=c("A", "C", "T", "G"))
  
  nt_df
  
}


nt_table_profile <- function(obj, main_suffix="all", ...) {
  
  if (inherits(obj, "GRanges"))
    obj <- prep_nt_table(obj, asfactor = T, ...)
  stopifnot(inherits(obj, "data.frame"))
  
  nt_freq_tab <- t(do.call(rbind, lapply(obj, table)))
  nt_freq_tab <- nt_freq_tab / colSums(nt_freq_tab)
  
  matplot(
    t(nt_freq_tab), type="l", ylab="freq", xlab="pos",
    main=paste0(
      "per-position nucleotide frequencies",
      "\nacross BS-centered 200bp regions (",main_suffix,")"),
    sub = paste(nrow(obj), "BSs"))
  abline(h=c(0,0.25,1), col="grey80")
  legend("topright", leg=rownames(nt_freq_tab),
         lty=1:5, col=1:6, text.col=1:6)
  
  return(invisible(nt_freq_tab))
  
}

  
nt_table_heatmap <- function(nt_df) {
  
  stopifnot(require(ggplot2) )
  stopifnot(require(reshape2))
  
  # convert rownames to factor, maintaining order seen in
  nt_df$ID <- factor(rownames(nt_df), levels=unique(rownames(nt_df)))
  
  # melt dataframe for use in geom_tile
  nt_df <- melt(nt_df, id.var = 'ID')
  colnames(nt_df) <- c("ID", "pos", "value")
  
  # ggplot
  ggplot(nt_df, aes_string(x="pos", y="ID")) + 
    geom_tile(aes_string(fill = "value"), colour = "grey50", width = 2, height = 2 )  +   
    scale_fill_manual(values = c(
      "A" = "darkgreen", "C" = "darkblue", "T" = "red", "G" = "yellow", "-" = "grey"))
}


avg_GC_outside_BS <- function(
    freq_tab, extension_width=100, avoid=1, FUN=median, scaling=100) {
  
  pos_idx <- c(
    1:(extension_width-avoid), 
    (ncol(freq_tab)-extension_width+1+avoid):ncol(freq_tab))
  
  scaling * FUN(colSums(freq_tab[c("C","G"),pos_idx]))
  
}

# EOF