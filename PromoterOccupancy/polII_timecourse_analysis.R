setwd("/work/gbioinfo/carlsara/grosshans/milou/polII_timecourse/")
library(QuasR)
library("BSgenome.Celegans.UCSC.ce10")
library(preprocessCore)

# initiate project
clObj <- makeCluster(20)
proj <- qAlign("samples.txt", "BSgenome.Celegans.UCSC.ce10", alignmentsDir="bam", clObj=clObj)

# Define TSS regions: +- 500b around TSS
exons <- read.delim("/work/gbioinfo/DB/WormBase/WS220/c_elegans.WS220.exons.tab", header=T, sep="\t", as.is=T)
exons <- exons[exons$chr!="chrM",] # remove mitochondrial genes
exons$chr <- as.factor(as.character(exons$chr)) # remove mitochondrial level

geneBodyStart <- aggregate(exons$start,list(exons$geneID),min)
geneBodyEnd <- aggregate(exons$end,list(exons$geneID),max)
geneBodyChr <- data.frame(do.call(rbind,strsplit(unique(paste(exons$geneID,exons$chr,sep=":")),":")))
geneBodyStrand <- data.frame(do.call(rbind,strsplit(unique(paste(exons$geneID,exons$strand,sep=":")),":")))
geneBodyStartEnd <- merge(geneBodyStart,geneBodyEnd,by.x=1,by.y=1)
geneBodyStartEndChr <- merge(geneBodyStartEnd,geneBodyChr,by.x=1,by.y=1)
geneBodyStartEndChrStrand <- merge(geneBodyStartEndChr,geneBodyStrand,by.x=1,by.y=1)
geneBody_gr <- GRanges(seqnames = Rle(geneBodyStartEndChrStrand[,4]), 
                       ranges = IRanges(geneBodyStartEndChrStrand[,2], end = geneBodyStartEndChrStrand[,3]), 
                       strand=geneBodyStartEndChrStrand[,5])
names(geneBody_gr) <- geneBodyStartEndChrStrand[,1]
geneBody_gr <- geneBody_gr + 500

tss <- geneBody_gr
start(tss) <- ifelse(strand(tss)=="+", start(tss), end(tss))
end(tss) <- start(tss)
tss_500 <- tss + 500

# Quantify
tss_500_counts <- qCount(proj, tss_500, clObj=clObj)

# Try quantile norm
edata = log2(tss_500_counts[,-1] + 8)
edata = edata[rowMeans(edata) > 3, ]
norm_edata = normalize.quantiles(as.matrix(edata[,13:24]))
# ==> final object: norm_edata

# (initially written by Sarah H.Carl; cut & slightly edited for brevity - AATS 13/09/2022)

# EOF
