writeLines(capture.output(date(),sessionInfo()), 'DataProcessing_1660R_sessioninfo.txt')


'
Process RNAseq read alignments

DESCRIPTION
-----------
Analysis was performed by: Milou Meeuse 

This code allows the user to count RNAseq reads and to quantify and plot gene expression

This script consists of the following steps: 
  (1) count reads
  (2) normalize read counts
  (3) plot gene expression
  (4) export expression files
  

EXPERIMENT DESCRIPTION
----------------------
Experiment was performed by: Milou Meeuse
Experiment identifier in GNomEx: 1660R
Sequencing protocol: TruSeq Illumina total RNA, HiSeq 50 cycle single end reads, HiSeq 2500, submitted 2018-06-29 

aim: investigate gene expression changes upon depletion of grh-1::degron::3xflag by adding 250 uM auxin (pl) or ethanol (min) at 21h of development at 20 degr in liquid. 
     RNA was sampled from 24.5h until 30h at 30 min intervals

REFERENCES
----------
https://bioconductor.org/packages/devel/bioc/vignettes/QuasR/inst/doc/QuasR.html


INPUT
-----
Data format: 
-sequencing reads: fastq.gz files on /work2/ggrossha/milou/RNAseq1660
-alignment proj: 2018-07-18 Proj.RData
-reference genome: BSgenome.Celegans.UCSC.ce10
-exon annotation: /work/gbioinfo/DB/WormBase/WS220/c_elegans.WS220.exons.tab


USAGE
-----
run qCount, normalize, plot and export



OUTPUT
------
-gene counts
-library size normalized, log2 transformed, pseudocount, expression: TLE 
-library size normalized, log2 transformed, pseudocount, mean normalized expression: TLED
-correlation plots





REQUIREMENTS
------------
This script requires the following packages to be installed:
-BSgenome_1.56.0
-NMF_0.22.0 
-BSgenome.Celegans.UCSC.ce10_1.4.0
-GenomicRanges_1.40.0
-IRanges_2.22.2


"'

setwd("/tungstenfs/scratch/ggrossha/meeumilo/scripts/ggrosshans_grh-1_paper")

library('BSgenome.Celegans.UCSC.ce10')
library('NMF')

#(1)---count reads---
clObj <- makeCluster(10)
exons <- read.delim('/work/gbioinfo/DB/WormBase/WS220/c_elegans.WS220.exons.tab',header=T,sep='\t',as.is=T)
head(exons)
exons <- exons[exons$chr!="chrM",] # remove mitochondrial genes
exons$chr <- as.factor(as.character(exons$chr)) # remove mitochondrial level
exons_gr <- GRanges(seqnames=exons$chr, ranges=IRanges(start=exons$start, end=exons$end), strand=exons$strand)
names(exons_gr) <- exons$geneID

gene_counts <- qCount(proj, exons_gr, orientation='same', clObj=clObj) #do both opposite and same strand and compare
gene_counts2 <- qCount(proj, exons_gr, orientation='opposite', clObj=clObj)
colSums(gene_counts)
colSums(gene_counts2)
write.table(gene_counts2,'Data/2020-07-17 gene_counts_1660R.tab',sep='\t',quote=FALSE)

#(2)---normalize read counts---
geneCounts <- read.table('Data/2020-07-17 gene_counts_1660R.tab',sep='\t')
gene_norm <- t(t(geneCounts[,-1])/colSums(geneCounts[,-1])*mean(colSums(geneCounts[,-1]))) #normalized to counts of each sample by scaling to the mean library size
TL <- log2(gene_norm + 8) #log2 transformed with pseudocount of 3
TLE <- TL[which(rowMeans(TL)>3),] #remove genes that are not expressed, i.e. TL=3
TLED <- TLE - rowMeans(TLE)

TL_pl <- TL[,c(seq(1,24,2))]
TL_min <- TL[,c(seq(2,24,2))]
TLE_pl <- TLE[,c(seq(1,24,2))]
TLE_min <- TLE[,c(seq(2,24,2))]
TLED_pl <- TLED[,c(seq(1,24,2))]
TLED_min <- TLED[,c(seq(2,24,2))]
TLED_min_pl <- cbind(TLED_min, TLED_pl) #normalization over both conditions allows to see differences
TLE_min_pl <- cbind(TLE_min, TLE_pl)


#(3)---plot gene expression---
aheatmap(cor(TLE_pl), Rowv=NA, Colv=NA) 
aheatmap(cor(TLE_min), Rowv=NA, Colv=NA) 
aheatmap(cor(TLE_pl, TLE_min), Rowv=NA, Colv=NA) 

plot(TLE_min[,1], TLE_pl[,1])
plot(TLE_min[,9], TLE_min[,10])
plot(TLE_min[,9], TLE_pl[,9])

pairs(TLE_pl, upper.panel = NULL) #sample 13p is off? maybe swapped with 12m?
aheatmap(cor(TLE[,c(1,3,5,7,9,11,13,15,17,18,21,23)]), Rowv = NA, Colv = NA) #swapped 13p with 12m
pairs(TLE_min, upper.panel = NULL) #sample 12m is off? maybe swapped with 13p?
aheatmap(cor(TLE[,c(2,4,6,8,10,12,14,16,19,20,22,24)]), Rowv = NA, Colv = NA) #swapped 12m with 113p
#note: sample 6m many unmapped reads


#(4)---export expression files---
write.csv(TLE_min_pl, file='Data/2020-07-22 TLEminpl_1660R.csv')
write.csv(TLE_min, file='Data/2020-07-22 TLEmin_1660R.csv')
write.csv(TLE_pl, file='Data/2020-07-22 TLEpl_1660R.csv')

write.csv(TLED_min_pl, file='Data/2020-07-22 TLEDminpl_1660R.csv')
write.csv(TLED_min, file='Data/2020-07-22 TLEDmin_1660R.csv')
write.csv(TLED_pl, file='Data/2020-07-22 TLEDpl_1660R.csv')

