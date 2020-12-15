#alignment was performed on 2018.11.27 QuasR version:	1.22.0. aligner:Rbowtie, version	1.22.0
#qCount was initially performed 2018.12.04, counted on wrong strand 
#qCount was performed again on 2019.02.06 

setwd('/work2/ggrossha/milou/RNAseq1779')

library('QuasR')
library('BSgenome.Celegans.UCSC.ce10')

##-align reads of RNA-seq to c.elegans genome ce10
clObj <- makeCluster(10)
proj <- qAlign('samples.txt', 'BSgenome.Celegans.UCSC.ce10', alignmentsDir = 'bam', clObj = clObj, splicedAlignment = TRUE)
exportwig <- qExportWig(proj, createBigWig = TRUE, binsize = 50)

#Quality Control report
qQCReport(proj, pdfFilename="QCReport.pdf")

##-create exon GRanges object: gene name, range of exon and strand
exons <- read.delim('/work/gbioinfo/DB/WormBase/WS220/c_elegans.WS220.exons.tab',header=T,sep='\t',as.is=T)
head(exons)
exons <- exons[exons$chr!="chrM",] # remove mitochondrial genes
exons$chr <- as.factor(as.character(exons$chr)) # remove mitochondrial level
exons_gr <- GRanges(seqnames=exons$chr, ranges=IRanges(start=exons$start, end=exons$end), strand=exons$strand)
names(exons_gr) <- exons$geneID

##-make qCount table: gene name, length gene and nr of counts
gene_counts <- qCount(proj, exons_gr, orientation='same', clObj=clObj) #do both opposite and same strand and compare
gene_counts2 <- qCount(proj, exons_gr, orientation='opposite', clObj=clObj)
plot(apply(gene_counts2[,-1],1,median), sd(gene_counts2[,-1])/rowMeans(gene_counts2[,-1]), xlim=c(0,500)) #median gene counts vs coeff of variation to exclude genes that are lowly expressed?
write.table(gene_counts2,'2019-02-06 gene_counts_2.tab',sep='\t',quote=FALSE)
