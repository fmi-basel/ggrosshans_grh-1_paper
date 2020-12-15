writeLines(capture.output(date(),sessionInfo()), 'Compare1660Rwith1779R_sessioninfo.txt')


'
Compare RNAseq 1660R with 1779R

DESCRIPTION
-----------
Analysis was performed by: Milou Meeuse 

This code allows the user to compare gene expression changes of RNAseq experiment 1660R with 1779R

This script consists of the following steps: 
  (1) load data
  (2) investigate similarity between time courses
  (3) plot gene expression changes of hits
  

EXPERIMENT DESCRIPTION
----------------------
Experiment was performed by: Milou Meeuse
Experiment identifier in GNomEx: 1660R and 1779R

aim: investigate gene expression changes upon depletion of grh-1::degron::3xflag by adding 250 uM auxin (pl) or ethanol (min) at 21h of development at 20 degr in liquid. 


REFERENCES
----------
https://bioconductor.org/packages/devel/bioc/vignettes/QuasR/inst/doc/QuasR.html
https://www.analyticsvidhya.com/blog/2016/03/pca-practical-guide-principal-component-analysis-python/


INPUT
-----
Data format: 
-TLE and TLED from 1660R (/tungstenfs/scratch/ggrossha/meeumilo/scripts/RNAseq_GRH-1degron/1660R/Data/) and 1779R (/tungstenfs/scratch/ggrossha/meeumilo/scripts/RNAseq_GRH-1degron/1779R/Data/)
-osc genes and their amplitude and phase (Meeuse et al, MSB 2020) (/tungstenfs/scratch/ggrossha/meeumilo/scripts/FullDevTC/2019-01-31 AllOsc_info_CIclass.csv)
-significant changing genes in 1779R (/tungstenfs/scratch/ggrossha/meeumilo/scripts/RNAseq_GRH-1degron/1779R/Data/2019-06-25 SignificantChangingGenes_RNAseq1779.csv)

USAGE
-----
Note that PCA showed that 28.5min has been swapped with 29pl during sample handling, which is part of this code.
To generate plots without exchanging timepoints, execute 1a and 1b, but not 1c
To generate plots with echanging timepoints, execute 1a, 1b, and 1c.


OUTPUT
------
- correlation heatmaps of 1660R and 1779R
- PCA of 1660R and 1779R
- gene expression heatmaps


REQUIREMENTS
------------
This script requires the following packages to be installed:
RColorBrewer_1.1-2  NMF_0.23.0          bigmemory_4.5.36    Biobase_2.48.0      BiocGenerics_0.34.0 cluster_2.1.0      
rngtools_1.5        pkgmaker_0.32.2     registry_0.5-1 


"'


setwd("/tungstenfs/scratch/ggrossha/meeumilo/scripts/ggrosshans_grh-1_paper")

library('NMF')

#(1)---load data---
#(a)1660R
TLEminpl_1660R <- as.matrix(read.csv('Data/2020-07-22 TLEminpl_1660R.csv', row.names = 1))
colnames(TLEminpl_1660R) <- c('24.5hmin', '25hmin', '25.5hmin', '26hmin','26.5hmin','27hmin','27.5hmin','28hmin','28.5hmin','29hmin','29.5hmin','30hmin','24.5hpl', '25hpl', '25.5hpl', '26hpl','26.5hpl','27hpl','27.5hpl','28hpl','28.5hpl','29hpl','29.5hpl','30hpl')
TLEDminpl_1660R <- as.matrix(read.csv('Data/2020-07-22 TLEDminpl_1660R.csv', row.names = 1))
colnames(TLEDminpl_1660R) <- c('24.5hmin', '25hmin', '25.5hmin', '26hmin','26.5hmin','27hmin','27.5hmin','28hmin','28.5hmin','29hmin','29.5hmin','30hmin','24.5hpl', '25hpl', '25.5hpl', '26hpl','26.5hpl','27hpl','27.5hpl','28hpl','28.5hpl','29hpl','29.5hpl','30hpl')
foldCh_1660R <- TLEminpl_1660R[,grep('pl',colnames(TLEminpl_1660R))]-TLEminpl_1660R[,grep('min',colnames(TLEminpl_1660R))] 
colnames(foldCh_1660R) <- sub('pl','',colnames(foldCh_1660R))

#(b)1779R
TLEminpl_1779R <- as.matrix(read.csv('Data/2019-06-25 TLEminpl_RNAseq1779.csv', row.names = 1)[,c(-11,-22)])
TLEDminpl_1779R<- as.matrix(read.csv('Data/2019-06-25 TLEDminpl_RNAseq1779.csv', row.names = 1)[,c(-11,-22)])
foldCh_1779R <- TLEminpl_1779R[,grep('pl',colnames(TLEminpl_1779R))]-TLEminpl_1779R[,grep('min',colnames(TLEminpl_1779R))] 
colnames(foldCh_1779R) <- sub('pl','',colnames(foldCh_1779R))

hits_1779R.info <- read.csv('Data/2020-10-26 DiffExpr_Annotated_RNAseq1779.csv', row.names = 1)
hits_1779R <- hits_1779R.info$WBid

osc_info <- read.csv('/tungstenfs/scratch/ggrossha/meeumilo/scripts/FullDevTC/2019-01-31 AllOsc_info_CIclass.csv', row.names = 1)

#(c)swap 28.5min with 29pl
TLEminpl_1660R <- TLEminpl_1660R[, c(1:8,22,10:21,9,23,24)] #swap 28.5min with 29pl
colnames(TLEminpl_1660R) <- c('24.5hmin', '25hmin', '25.5hmin', '26hmin','26.5hmin','27hmin','27.5hmin','28hmin','28.5hmin','29hmin','29.5hmin','30hmin','24.5hpl', '25hpl', '25.5hpl', '26hpl','26.5hpl','27hpl','27.5hpl','28hpl','28.5hpl','29hpl','29.5hpl','30hpl')
TLEDminpl_1660R <- TLEDminpl_1660R[, c(1:8,22,10:21,9,23,24)] #swap 28.5min with 29pl
colnames(TLEDminpl_1660R) <- c('24.5hmin', '25hmin', '25.5hmin', '26hmin','26.5hmin','27hmin','27.5hmin','28hmin','28.5hmin','29hmin','29.5hmin','30hmin','24.5hpl', '25hpl', '25.5hpl', '26hpl','26.5hpl','27hpl','27.5hpl','28hpl','28.5hpl','29hpl','29.5hpl','30hpl')
foldCh_1660R <- TLEminpl_1660R[,13:24]-TLEminpl_1660R[,1:12] 
colnames(foldCh_1660R) <- c('24.5h','25h','25.5h','26h','26.5h','27h','27.5h','28h','28.5h','29h','29.5h','30h')

#write.csv(TLEminpl_1660R, '2020-10-27 Compare1660with1779/2020-10-27 TLEminpl_1660R_TPswapped.csv')
#write.csv(TLEDminpl_1660R, '2020-10-27 Compare1660with1779/2020-10-27 TLEDminpl_1660R_TPswapped.csv')
#write.csv(foldCh_1660R, '2020-10-27 Compare1660with1779/2020-10-27 foldChange_1660R_TPswapped.csv')


#(2)---similarity between TCs---

#combine TCs in one matrix
genes_overlap <- rownames(TLEminpl_1660R)[rownames(TLEminpl_1660R) %in% rownames(TLEminpl_1779R)]
TLEminpl_comb <- cbind(TLEminpl_1779R[genes_overlap,], TLEminpl_1660R[genes_overlap,])
TLEmin_comb <- cbind(TLEminpl_1779R[genes_overlap,grep('min',colnames(TLEminpl_1779R))], TLEminpl_1660R[genes_overlap,grep('min',colnames(TLEminpl_1660R))])
TLEpl_comb <- cbind(TLEminpl_1779R[genes_overlap,grep('pl',colnames(TLEminpl_1779R))], TLEminpl_1660R[genes_overlap,grep('pl',colnames(TLEminpl_1660R))])
TLEDminpl_comb <- cbind(TLEDminpl_1779R[genes_overlap,], TLEDminpl_1660R[genes_overlap,])
TLEDmin_comb <- cbind(TLEDminpl_1779R[genes_overlap,grep('min',colnames(TLEDminpl_1779R))], TLEDminpl_1660R[genes_overlap,grep('min',colnames(TLEDminpl_1660R))])
TLEDpl_comb <- cbind(TLEDminpl_1779R[genes_overlap,grep('pl',colnames(TLEDminpl_1779R))], TLEDminpl_1660R[genes_overlap,grep('pl',colnames(TLEDminpl_1660R))])

#(a)-correlation of all genes and osc genes-
aheatmap(cor(TLEmin_comb), Rowv = NA, Colv = NA) #28.5min (1660R) is off, timing of TCs is not equal
aheatmap(cor(TLEpl_comb), Rowv = NA, Colv = NA) #29pl (1660R) is off, timing of TCs is not equal
aheatmap(cor(TLEminpl_comb[rownames(osc_info),]), Rowv = NA, Colv = NA) #28.5min and 29pl swapped?
aheatmap(cor(TLEmin_comb[rownames(osc_info),]), Rowv = NA, Colv = NA)
aheatmap(cor(TLEpl_comb[rownames(osc_info),]), Rowv = NA, Colv = NA)

#(b)-PCA on TLE of osc genes-
  #note: because PCA is on TLE, PC1 represents the mean expression
pca_1660R <- prcomp(TLEminpl_1660R[rownames(osc_info),])
pca_1660R_var <- pca_1660R$sdev^2
barplot(100*pca_1660R_var/sum(pca_1660R_var))

pca_1779R <- prcomp(TLEminpl_1779R[rownames(osc_info),])
pca_1779R_var <- pca_1779R$sdev^2
barplot(100*pca_1779R_var/sum(pca_1779R_var))

pca_comb <- prcomp(TLEminpl_comb[rownames(osc_info),])
pca_comb_var <- pca_comb$sdev^2
barplot(100*pca_comb_var/sum(pca_comb_var))

par(mfrow=c(1,2))
plot(pca_1660R$rotation[,3], pca_1660R$rotation[,2], col = c(rep('purple',12),rep('orange1',12)), main='1660R', xlab='PC3', ylab='PC2');
lines(pca_1660R$rotation[1:12,3], pca_1660R$rotation[1:12,2], col='purple'); 
lines(pca_1660R$rotation[13:24,3], pca_1660R$rotation[13:24,2], col='orange1'); 
text(pca_1660R$rotation[1:12,3], pca_1660R$rotation[1:12,2], colnames(TLEDminpl_1660R)[1:12], col='purple', pos = 4, cex=0.75);
text(pca_1660R$rotation[13:24,3], pca_1660R$rotation[13:24,2], colnames(TLEDminpl_1660R)[13:24], col='orange1', pos = 4, cex=0.75);

plot(pca_1779R$rotation[,3], pca_1779R$rotation[,2], col = c(rep('orangered1',10),rep('blue',10)), main='1779R', xlab='PC3', ylab='PC2');
lines(pca_1779R$rotation[1:10,3], pca_1779R$rotation[1:10,2], col='orangered1'); 
lines(pca_1779R$rotation[11:20,3], pca_1779R$rotation[11:20,2], col='blue'); 
text(pca_1779R$rotation[1:10,3], pca_1779R$rotation[1:10,2], colnames(TLEDminpl_1779R)[1:10], col='orangered1', pos = 4, cex=0.75);
text(pca_1779R$rotation[11:20,3], pca_1779R$rotation[11:20,2], colnames(TLEDminpl_1779R)[11:22], col='blue', pos = 4, cex=0.75);

par(mfrow=c(1,1))
plot(pca_comb$rotation[,3], pca_comb$rotation[,2], col = c(rep('orangered1',10),rep('blue',10),rep('purple',12), rep('orange1',12)),main='1660R and 1779R',xlab='PC3', ylab='PC2'); 
lines(pca_comb$rotation[21:32,3], pca_comb$rotation[21:32,2], col='purple'); 
lines(pca_comb$rotation[33:44,3], pca_comb$rotation[33:44,2], col='orange1'); 
text(pca_comb$rotation[21:32,3], pca_comb$rotation[21:32,2], colnames(TLEDminpl_comb)[21:32], col='purple', pos = 4, cex=0.75);
text(pca_comb$rotation[33:44,3], pca_comb$rotation[33:44,2], colnames(TLEDminpl_comb)[33:44], col='orange1', pos = 4, cex=0.75);
lines(pca_comb$rotation[1:10,3], pca_comb$rotation[1:10,2], col='orangered1'); 
lines(pca_comb$rotation[11:20,3], pca_comb$rotation[11:20,2], col='blue'); 
text(pca_comb$rotation[1:10,3], pca_comb$rotation[1:10,2], colnames(TLEDminpl_comb)[1:10], col='orangered1', pos = 4, cex=0.75);
text(pca_comb$rotation[11:20,3], pca_comb$rotation[11:20,2], colnames(TLEDminpl_comb)[11:20], col='blue', pos = 4, cex=0.75);


#(3)---plot gene expression changes of hits---
h1 <- hclust(dist(TLEDminpl_1779R[hits_1779R,]))
hits_cl <- hits_1779R[h1$order] #sort by hierarchical clustering
par(mfrow=c(1,4))
breaks <- c(seq(-5,5, length=256))
aheatmap(pmax(pmin(TLEDminpl_1779R[hits_cl,],5),-5), Colv = NA, Rowv = NA, col="-RdBu:256", breaks=breaks, main ='1779R')
aheatmap(pmax(pmin(TLEDminpl_1660R[hits_cl,],5),-5), Colv = NA, Rowv = NA, col="-RdBu:256", breaks=breaks, main='1660R')
aheatmap(pmax(pmin(foldCh_1779R[hits_cl,],5),-5), Colv = NA, Rowv = NA, col="-RdBu:256", breaks=breaks, main ='1779R')
aheatmap(pmax(pmin(foldCh_1660R[hits_cl,],5),-5), Colv = NA, Rowv = NA, col="-RdBu:256", breaks=breaks, main ='1660R')

h2 <- hclust(dist(foldCh_1779R[hits_1779R,]))
hits_cl2 <- hits_1779R[h2$order] #sort by hierarchical clustering
par(mfrow=c(1,4))
breaks <- c(seq(-5,5, length=256))
aheatmap(pmax(pmin(foldCh_1779R[hits_cl2,],5),-5), Colv = NA, Rowv = NA, col="-RdBu:256", breaks=breaks, main ='1779R')
aheatmap(pmax(pmin(foldCh_1660R[hits_cl2,],5),-5), Colv = NA, Rowv = NA, col="-RdBu:256", breaks=breaks, main ='1660R')
aheatmap(pmax(pmin(TLEDminpl_1779R[hits_cl2,],5),-5), Colv = NA, Rowv = NA, col="-RdBu:256", breaks=breaks, main ='1779R')
aheatmap(pmax(pmin(TLEDminpl_1660R[hits_cl2,],5),-5), Colv = NA, Rowv = NA, col="-RdBu:256", breaks=breaks, main='1660R')

nonHits<-rownames(foldCh_1660R)[!(rownames(foldCh_1660R) %in% hits_1779R)]
h3 <- hclust(dist(foldCh_1660R[nonHits,]))
hits_cl3 <- nonHits[h3$order] #sort by hierarchical clustering
aheatmap(pmax(pmin(foldCh_1660R[hits_cl3,],5),-5), Colv = NA,  Rowv = NA, col="-RdBu:256", breaks=breaks, main ='1660R - non-hits')
aheatmap(pmax(pmin(TLEDminpl_1660R[hits_cl3,],5),-5), Colv = NA,  Rowv = NA, col="-RdBu:256", breaks=breaks, main ='1660R - non-hits')

nonHits1779R<-rownames(foldCh_1779R)[!(rownames(foldCh_1779R) %in% hits_1779R)]
h4 <- hclust(dist(foldCh_1779R[nonHits1779R,]))
hits_cl4 <- nonHits1779R[h4$order] #sort by hierarchical clustering
aheatmap(pmax(pmin(foldCh_1779R[hits_cl4,],5),-5), Colv = NA,  Rowv = NA, col="-RdBu:256", breaks=breaks, main ='1779R - non-hits')
aheatmap(pmax(pmin(TLEDminpl_1779R[hits_cl4,],5),-5), Colv = NA,  Rowv = NA, col="-RdBu:256", breaks=breaks, main ='1779R - non-hits')










