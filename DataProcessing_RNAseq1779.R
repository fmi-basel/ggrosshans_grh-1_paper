writeLines(capture.output(sessionInfo()), "2019-06-25 DataProcessing_RNAseq1779_sessionInfo.txt")

setwd('/tungstenfs/scratch/ggrossha/meeumilo/scripts/RNAseq_GRH-1degron/1779R')

library('NMF')

##-normalize gene counts
geneCounts <- read.table('Data/2019-02-06 gene_counts_2.tab',sep='\t')
gene_norm <- t(t(geneCounts[,-1])/colSums(geneCounts[,-1])*mean(colSums(geneCounts[,-1]))) #normalized to counts os each sample by scaling to the mean library size
TL <- log2(gene_norm + 8) #log2 transformed with pseudocount of 3
TLE <- TL[which(rowMeans(TL)>3),] #remove genes that are not expressed, i.e. TL=3
TLED <- TLE - rowMeans(TLE)

##-split in + and - auxin condition
TL_pl <- TL[,c(1,3,5,7,9,11,13,15,17,19,21)]
TL_min <- TL[,c(2,4,6,8,10,12,14,16,18,20,22)]
TLE_pl <- TLE[,c(1,3,5,7,9,11,13,15,17,19,21)]
TLE_min <- TLE[,c(2,4,6,8,10,12,14,16,18,20,22)]
TLED_pl <- TLED[,c(1,3,5,7,9,11,13,15,17,19,21)]
TLED_min <- TLED[,c(2,4,6,8,10,12,14,16,18,20,22)]
TLED_min_pl <- cbind(TLED_min, TLED_pl) #normalization over both conditions allows to see differences
TLE_min_pl <- cbind(TLE_min, TLE_pl)

write.csv(TLE_min_pl, file='Data/2019-06-25 TLEminpl_RNAseq1779.csv')
write.csv(TLE_min, file='Data/2019-06-25 TLEmin_RNAseq1779.csv')
write.csv(TLE_pl, file='Data/2019-06-25 TLEpl_RNAseq1779.csv')

write.csv(TLED_min_pl, file='Data/2019-06-25 TLEDminpl_RNAseq1779.csv')
write.csv(TLED_min, file='Data/2019-06-25 TLEDmin_RNAseq1779.csv')
write.csv(TLED_pl, file='Data/2019-06-25 TLEDpl_RNAseq1779.csv')


##-plot normalized log2 transformed gene counts of one sample vs other and create correlation heatmap with pseudocount 3
par(mfrow=c(1,3))
aheatmap(cor(TL_pl), Rowv=NA, Colv=NA) 
aheatmap(cor(TL_min), Rowv=NA, Colv=NA) 
aheatmap(cor(TL_pl, TL_min), Rowv=NA, Colv=NA) 

plot(TL[,1],TL[,2],ylab='22h_pl',xlab='22h_min') 
plot(TL[,1],TL[,3],ylab='22h_pl',xlab='23h_pl') 
plot(TL[,1],TL[,21],ylab='22h_pl',xlab='32h_pl')

par(mfrow=c(2,3))
plot(TLE_min[,1], TLE_pl[,1], ylim=c(3,21), xlim=c(3,21), main='22h', xlab='-auxin', ylab='+auxin')
plot(TLE_min[,2], TLE_pl[,2], ylim=c(3,21), xlim=c(3,21), main='23h', xlab='-auxin', ylab='+auxin')
plot(TLE_min[,3], TLE_pl[,3], ylim=c(3,21), xlim=c(3,21), main='24h', xlab='-auxin', ylab='+auxin')
plot(TLE_min[,4], TLE_pl[,4], ylim=c(3,21), xlim=c(3,21), main='25h', xlab='-auxin', ylab='+auxin')
plot(TLE_min[,5], TLE_pl[,5], ylim=c(3,21), xlim=c(3,21), main='26h', xlab='-auxin', ylab='+auxin')
plot(TLE_min[,6], TLE_pl[,6], ylim=c(3,21), xlim=c(3,21), main='27h', xlab='-auxin', ylab='+auxin')

plot(TLE_min[,7], TLE_pl[,7], ylim=c(3,21), xlim=c(3,21), main='28h', xlab='-auxin', ylab='+auxin')
plot(TLE_min[,8], TLE_pl[,8], ylim=c(3,21), xlim=c(3,21), main='29h', xlab='-auxin', ylab='+auxin')
plot(TLE_min[,9], TLE_pl[,9], ylim=c(3,21), xlim=c(3,21), main='30h', xlab='-auxin', ylab='+auxin')
plot(TLE_min[,10], TLE_pl[,10], ylim=c(3,21), xlim=c(3,21), main='31h', xlab='-auxin', ylab='+auxin')
plot(TLE_min[,11], TLE_pl[,11], ylim=c(3,21), xlim=c(3,21), main='32h', xlab='-auxin', ylab='+auxin')

