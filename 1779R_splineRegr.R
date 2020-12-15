writeLines(capture.output(date(),sessionInfo()), '1779R_splineRegr_sessioninfo.txt')


'
Differential expression analysis using spline regression

DESCRIPTION
-----------
Analysis was performed by: Milou Meeuse 

This code allows the user to apply spline regression on delta expression between conditions to identify differentially expressed genes without having replicate time courses

This script consists of the following steps: 
  (1) load data
  (2) differential gene expression using regression with bspline
  (3) select differentially expressed genes
  (4) phase specificity
  (5) export
  

EXPERIMENT DESCRIPTION
----------------------
Experiment was performed by: Milou Meeuse
Experiment identifier in GNomEx: 1779R

aim: investigate gene expression changes upon depletion of grh-1::degron::3xflag by adding 250 uM auxin (pl) or ethanol (min) at 21h of development at 20 degr in liquid. 


REFERENCES
----------


INPUT
-----
- gene expression in pl and min condition: Data/2019-06-25 TLEminpl_RNAseq1779.csv
- oscillating genes according to Meeuse et al, 2020: /tungstenfs/scratch/ggrossha/meeumilo/scripts/FullDevTC/2019-01-31 AllOsc_info_CIclass.csv
- gene annotation ce10: /tungstenfs/scratch/ggrossha/meeumilo/scripts/geneAnnotation/Celegans_WS220_WBid-geneName-transcrName-type.csv


OUTPUT
------
- list of differentially expressed genes
- heatmaps with differential gene expression
- histograms of phase specificity


REQUIREMENTS
------------
RColorBrewer_1.1-2  NMF_0.23.0          bigmemory_4.5.36    Biobase_2.48.0      BiocGenerics_0.34.0 cluster_2.1.0      
rngtools_1.5        pkgmaker_0.31.1     registry_0.5-1   


"'

setwd("/tungstenfs/scratch/ggrossha/meeumilo/scripts/ggrosshans_grh-1_paper")
library('splines')
library(NMF)


#(1)---load data---
TLE_min_pl <- as.matrix(read.csv('Data/2019-06-25 TLEminpl_RNAseq1779.csv', row.names = 1))[,c(1:10,12:21)] #omit TP11, as development is severely affected, see PCA in 1779R_AlignmentOfConditions.docx
TLE_min <- TLE_min_pl[,1:10] 
TLE_pl <- TLE_min_pl[,11:20] 
TLED_min_pl <- TLE_min_pl - rowMeans(TLE_min_pl)
FullDevTC_AmplPhase <- read.csv('/tungstenfs/scratch/ggrossha/meeumilo/scripts/FullDevTC/2019-01-31 AllOsc_info_CIclass.csv', header=T,sep=',', row.names = 1)
FullDevTC_AmplPhase.sort <- FullDevTC_AmplPhase[sort.int(FullDevTC_AmplPhase$Phase, index.return = T)$ix,]
FullDevTC_Osc <- rownames(FullDevTC_AmplPhase)


#(2)---differential gene expression using regression with bspline---
delta <- TLE_pl - TLE_min
bsplineb <- bs(1:ncol(delta),df=5) # generate B-spline basis matrix for cubic spline
matplot(bsplineb, type='l') 

#plot example gene
par(mfrow=c(2,1))
plot(seq(22,31,by=1),TLE_min['WBGene00014009',], type='l', ylim=c(4,15), ylab='lips-6 expression (log2)', xlab= 'Time of development (h)');lines(seq(22,31,by=1),TLE_pl['WBGene00014009',], col='turquoise3')
plot(seq(22,31,by=1),delta['WBGene00014009',], col='orangered', type='l', ylab='expression fold change', xlab= 'Time of development (h)') #example gene

#fit on example gene
delta.sel <- delta['WBGene00001707',]
fit.sel <- lm(delta.sel ~ bsplineb)
plot(delta.sel);lines(fit.sel$fitted.values)
matplot(cbind(fit.sel$coefficients[2]*bsplineb[,1], fit.sel$coefficients[3]*bsplineb[,2],fit.sel$coefficients[4]*bsplineb[,3],fit.sel$coefficients[5]*bsplineb[,4]), type='l')
plot(fit.sel$coefficients[2]*bsplineb[,1] + fit.sel$coefficients[3]*bsplineb[,2] + fit.sel$coefficients[4]*bsplineb[,3] + fit.sel$coefficients[5]*bsplineb[,4], type='l')

#all genes
fit <- lm(t(delta)~bsplineb) 
fit.sum <- summary(fit)
fit.fstat <- t(sapply(fit.sum,function(x){x$fstatistic}))
fit.pval <- pf(fit.fstat[,1], fit.fstat[,2], fit.fstat[,3], lower.tail = F)
fit.pvalAdj <- p.adjust(fit.pval, method = 'fdr')
names(fit.pvalAdj) <- sub('Response ','',names(fit.pvalAdj))


#(3)---select differentially expressed genes---
plot(-log10(fit.pval), -log10(fit.pvalAdj)); abline(0,1)
hist(-log10(fit.pvalAdj))

delta.fit <- t(fit$fitted.values)
delta.fit.Abs <- abs(delta.fit)
plot(-log10(fit.pvalAdj),rowSums(delta.fit.Abs))

hits_spline <- names(which(fit.pvalAdj<0.05 & rowSums(delta.fit.Abs>0.5)>=3))
length(hits_spline)
sum(hits_spline %in% FullDevTC_Osc)
sum(hits_spline %in% FullDevTC_Osc)/length(hits_spline)

hits_spline.osc <- hits_spline[hits_spline %in% FullDevTC_Osc]
hits_spline.osc.sort <- rownames(FullDevTC_AmplPhase.sort[rownames(FullDevTC_AmplPhase.sort) %in% hits_spline.osc,])

h1<-hclust(dist(delta[hits_spline,],method="euclidean"), method="complete")
par(mfrow=c(1,3))
breaks <- seq(-5,5,length.out = 256)
aheatmap(pmax(pmin(delta[hits_spline[h1$order],],5),-5), Rowv = NA, Colv = NA, col="-RdBu:256",annRow = list(osc=factor(hits_spline %in% FullDevTC_Osc)), annColors = list(osc=c('olivedrab1','purple4')), labRow = NA, breaks=breaks, labCol = c(seq(22,31,by=1)), main='gene expression fold change (n=633)')
aheatmap(pmax(pmin(TLED_min_pl[hits_spline[h1$order],1:10],5),-5), Rowv = NA, Colv = NA, col="-RdBu:256", breaks=breaks, labRow = NA, labCol = c(seq(22,31,by=1)), main='normalized gene expression (n=633) - vehicle')
aheatmap(pmax(pmin(TLED_min_pl[hits_spline[h1$order],11:20],5),-5), Rowv = NA, Colv = NA, col="-RdBu:256", breaks=breaks, labRow = NA, labCol = c(seq(22,31,by=1)), main='- control')

aheatmap(pmax(pmin(delta[hits_spline.osc.sort,],5),-5),Rowv = NA,  Colv = NA, col="-RdBu:256", breaks=breaks, labRow = NA, labCol = c(seq(22,31,by=1)), main='gene expression fold change (n=542)')
aheatmap(pmax(pmin(TLED_min_pl[hits_spline.osc.sort,1:10],5),-5), Rowv = NA, Colv = NA, col="-RdBu:256", breaks=breaks, labRow = NA, labCol = c(seq(22,31,by=1)), main='normalized gene expression (n=542) - vehicle')
aheatmap(pmax(pmin(TLED_min_pl[hits_spline.osc.sort,11:20],5),-5), Rowv = NA, Colv = NA, col="-RdBu:256", breaks=breaks, labRow = NA, labCol = c(seq(22,31,by=1)), main='- control')

par(mfrow=c(8,9),mar=c(1.8, 1.8, 0.5, 0.5))
for (i in 1:72){
  gene <- hits_spline[i]
  plot(TLE_pl[gene,], col='red', ylab = 'expr', type='l',ylim=c(min(c(TLE_pl[gene,],TLE_min[gene,])), max(c(TLE_pl[gene,],TLE_min[gene,])))); lines(TLE_min[gene,], col='black')
}

for (i in 1:72){
  gene <- hits_spline[i]
  plot(delta[gene,], ylab = 'delta expr', ylim=c(min(c(delta[gene,],delta.fit[gene,])), max(c(delta[gene,],delta.fit[gene,])))); lines(delta.fit[gene,], col='red')
}


#(4)---phase specificity---
par(mfrow=c(1,2))
hist.phase.osc <- hist(FullDevTC_AmplPhase[,2], breaks=30, freq = F, main='oscillating genes')
hist.phase.hitsOsc <- hist(FullDevTC_AmplPhase[hits_spline.osc,2], breaks=30, freq = F, main='oscillating hits') 

hits_spline.osc.up <- names(which(rowSums(delta[hits_spline.osc,])>=0))
hits_spline.osc.down <- names(which(rowSums(delta[hits_spline.osc,])<0))
hist.phase.hitsOsc.up <- hist(FullDevTC_AmplPhase[hits_spline.osc.up,2], breaks=30,freq = F, main='up oscillating hits') 
hist.phase.hitsOsc.down <- hist(FullDevTC_AmplPhase[hits_spline.osc.down,2], breaks=30,freq = F, main='down oscillating hits') 

#fraction
frac.phase.osc <- hist.phase.osc$counts/nrow(FullDevTC_AmplPhase) #fraction of osc genes in each bin (all osc genes)
pred.phase.hitsOsc.up <- frac.phase.osc*length(hits_spline.osc.up)
pred.phase.hitsOsc.down <- frac.phase.osc*length(hits_spline.osc.down)
fc.phase.hitsOsc.up <- (hist.phase.hitsOsc.up$counts+8)/(pred.phase.hitsOsc.up+8) #	To	minimize	the	large	enrichments	that	would	otherwise	be	caused	by	tissues	with	small	number	of	genes,	we	added	a	pseudocount	of	8 before	calculating	the	actual	ratio
fc.phase.hitsOsc.down <- (hist.phase.hitsOsc.down$counts+8)/(pred.phase.hitsOsc.down+8) #	To	minimize	the	large	enrichments	that	would	otherwise	be	caused	by	tissues	with	small	number	of	genes,	we	added	a	pseudocount	of	8 before	calculating	the	actual	ratio

#to get the directionallity in the plot, adjust again with manually setting the axis
barplot(fc.phase.hitsOsc.up-1, yaxt='n', ylim=c(-0.4,1), ylab='fold change over all oscillating genes', col='tomato3', names.arg=hist.phase.hitsOsc.up$mids, las=2, main='Upregulated (n=243)')
axis(2, at=seq(-0.4,1, by=0.2), labels=seq(-0.4,1, by=0.2)+1) 
barplot(fc.phase.hitsOsc.down-1, yaxt='n', ylim=c(-0.4,1), ylab='fold change over all oscillating genes', col='navyblue', names.arg=hist.phase.hitsOsc.down$mids, las=2, main='Downregulated (n=299)')
axis(2, at=seq(-0.4,1, by=0.2), labels=seq(-0.4,1, by=0.2)+1)

#grh-1 at 269.2 degrees
#molt at 


#(5)---export---
geneAnnot <- read.csv('/tungstenfs/scratch/ggrossha/meeumilo/scripts/geneAnnotation/Celegans_WS220_WBid-geneName-transcrName-type.csv', row.names=1)
rownames(geneAnnot) <- geneAnnot$WBID
geneAnnot <- geneAnnot[,-1]
hits_spline.annot <- geneAnnot[hits_spline,]
hits_spline.annot$osc <-  hits_spline %in% FullDevTC_Osc
hits_spline.annot <- merge(hits_spline.annot,FullDevTC_AmplPhase[rownames(FullDevTC_AmplPhase) %in% hits_spline,], all = T,by='row.names')
colnames(hits_spline.annot)[1] <- 'WBid'
hits_spline.delta <- ifelse(rowSums(delta[hits_spline,])>=0, 'up','down')
hits_spline.annot$delta <- hits_spline.delta[hits_spline.annot$WBid]
write.csv(hits_spline.annot, 'Data/2020-10-26 DiffExpr_Annotated_RNAseq1779.csv')



