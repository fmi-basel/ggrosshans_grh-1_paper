# Analyse mRNA-seq TC data
# (initially written by Milou WM Meeuse; cut & slightly edited for brevity by AATS 13/09/2022)

# set up
setwd('/work2/ggrossha/milou/RNAseq1314')
library('QuasR')
library('BSgenome.Celegans.UCSC.ce10')


# align reads of RNA-seq to c.elegans genome
clObj <- makeCluster(10)
proj <- qAlign('samples.txt', 'BSgenome.Celegans.UCSC.ce10', alignmentsDir = 'bam', clObj = clObj, splicedAlignment = TRUE)


# create exon GRanges object: gene name, range of exon and strand
exons <- read.delim('/work/gbioinfo/DB/WormBase/WS220/c_elegans.WS220.exons.tab',header=T,sep='\t',as.is=T)
exons <- exons[exons$chr!="chrM",] # remove mitochondrial genes
exons$chr <- as.factor(as.character(exons$chr)) # remove mitochondrial level
exons_gr <- GRanges(seqnames=exons$chr, ranges=IRanges(start=exons$start, end=exons$end), strand=exons$strand)
names(exons_gr) <- exons$geneID


# make qCount table: gene name, length gene and nr of counts
gene_counts <- qCount(proj, exons_gr, orientation='same', clObj=clObj) # do gene counts for both orientations to make sure you choose the right one


# normalize counts for library size and log2 transform
gene_norm <- t(t(gene_counts[,-1])/colSums(gene_counts[,-1])*mean(colSums(gene_counts[,-1])))
TL <- log2(gene_norm + 8)
# ==> Final object exported for firther use: TL




# EOF
