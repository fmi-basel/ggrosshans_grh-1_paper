writeLines(capture.output(date(),sessionInfo()), "Matrisome_Tissues_sessionInfo.txt")


'
Investigate tissue, cell type and matrisome enrichment of grh-1 degron hits 

DESCRIPTION
-----------
This code allows the user to investigate tissue enrichment (Cao), cell type enrichment (Cao) and matrisome enrichment (Teuscher) of a list of genes

This script consists of the following steps: 
  (1) max tissue, i.e. tissue with the highest expression and the ratio with the second tissue >=5 and qval<0.05
  (2) single cell expression: expression of selected genes in single cell data for each cell type
  (3) matrisome: classify selected genes according to matrisome categories
  

REFERENCES
----------
tissue and single cell:
Cao et al, Science, 2017, https://science.sciencemag.org/content/357/6352/661.long
Table S4: https://science.sciencemag.org/content/suppl/2017/08/16/357.6352.661.DC1?_ga=2.176460325.1014943348.1589951240-1324506074.1556476798

matrisome:
Teuscher, 2019: https://www.sciencedirect.com/science/article/pii/S2590028518300012
http://ec2-3-120-159-30.eu-central-1.compute.amazonaws.com:3838/ubuntu/ecm_analyzer/ (accessed on 26.10.2020)

fold change: calculated as example on p12: http://pantherdb.org/help/PANTHER_user_manual.pdf
predicted = nr of genes in that class/total nr of genes in the reference list (genome) * nr of genes of interest
enrichment =  genes of interest in class + 8 / predicted + 8 (To	minimize	the	large	enrichments	that	would	otherwise	be	caused	by	tissues	with	small	number	of	genes,	we	added	a	pseudocount	of	8 before	calculating	the	actual	ratio) 

INPUT
-----
Data format: 
- Cao-tissues: list of tissue (max tissue, 1st tissue) with max expression, 2nd tissue with expression, ratio of 1st and 2nd tissue and pval/qvalues
- Cao-single cell: matrix of transcripts per million of each gene in each cell type identified by scRNAseq
- genes: list of genes to investigate
- matrisome results: for each gene in the gene list provided the division and category in the matrisome


OUTPUT
------
- enrichment for each: tissues, cell type, matrisome in a barplot


REQUIREMENTS
------------
RColorBrewer_1.1-2  NMF_0.23.0          bigmemory_4.5.36    Biobase_2.48.0      BiocGenerics_0.34.0 cluster_2.1.0      
rngtools_1.5        pkgmaker_0.32.2     registry_0.5-1  



"'


setwd("/tungstenfs/scratch/ggrossha/meeumilo/scripts/ggrosshans_grh-1_paper")
library(NMF)

#(1)---Max tissue---
Cao.tis <- read.csv('/tungstenfs/scratch/ggrossha/meeumilo/scripts/TissueExpr/2017_Cao_Shendure_single_cell_Tissue_information_genes_full.csv', row.names = 1)
Cao.tisS <- Cao.tis[which(Cao.tis$ratio>=5 & Cao.tis$qval<0.05),]

genesAnnot <- read.csv('Data/2020-10-26 DiffExpr_Annotated_RNAseq1779.csv', row.names = 1)
row.names(genesAnnot) <- genesAnnot$WBid
genes <- genesAnnot$WBid
sum(rownames(Cao.tisS) %in% genes)

max.tis <- data.frame(allGenes=as.matrix(table(Cao.tisS$max.tissue)), hits=as.matrix(table(Cao.tisS[rownames(Cao.tisS) %in% genes,'max.tissue'])))
max.tis$pred <- (max.tis$allGenes/20000)*length(genes) #nr predicted in hits based on all genes
max.tis$FC <- (max.tis$hits+8)/(max.tis$pred+8) 
max.tis$FC_cor <- max.tis$FC-1 #to get the directionallity in the plot, adjust again with manually setting the axis
barplot(t(max.tis$FC_cor), beside=TRUE, horiz = TRUE, las=TRUE, xlim=c(-1,5),xaxt = "n", xlab='Fold enrichment', names.arg = rownames(max.tis), main=paste('n=',sum(rownames(Cao.tisS) %in% genes), sep=''))
axis(1, at=seq(-1,5, by=1), labels=seq(-1,5, by=1)+1)


#(2)---single cell expression---
Cao.sc <- read.csv('/tungstenfs/scratch/ggrossha/meeumilo/scripts/TissueExpr/Cao2017_tableS4.csv', row.names=1)
Cao.scL <- log2(as.matrix(Cao.sc[,-1])+8)

Cao.scLS <- Cao.scL[rownames(Cao.scL) %in% genes,]
genesAnnot.Cao.sc <- genesAnnot[rownames(Cao.scLS),]

#plot all
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
aheatmap(t(Cao.scLS), col=jet.colors(256),Rowv = FALSE, Colv = FALSE,annCol = list(osc=factor(genesAnnot.Cao.sc$osc), delta=factor(genesAnnot.Cao.sc$delta)), annColors = list(osc=c('olivedrab1','purple4'), delta=c('navyblue','tomato3'),grh1=c('grey90', 'black')), annRow = list(grh1=Cao.scLS['WBGene00001707',]))

#plot up and down separately
par(mfrow=c(1,2))
aheatmap(t(Cao.scLS[genesAnnot.Cao.sc$delta=='up',]), col=jet.colors(256),Rowv = NA, Colv = FALSE,annCol = list(osc=factor(genesAnnot.Cao.sc$osc[genesAnnot.Cao.sc$delta=='up']), delta=factor(genesAnnot.Cao.sc$delta[genesAnnot.Cao.sc$delta=='up'])), annColors = list(osc=c('olivedrab1','purple4'), delta=c('tomato3','navyblue'),grh1=c('grey90', 'black')), annRow = list(grh1=Cao.scLS['WBGene00001707',]))
aheatmap(t(Cao.scLS[genesAnnot.Cao.sc$delta=='down',]), col=jet.colors(256),Rowv = NA, Colv = FALSE,annCol = list(osc=factor(genesAnnot.Cao.sc$osc[genesAnnot.Cao.sc$delta=='down']), delta=factor(genesAnnot.Cao.sc$delta[genesAnnot.Cao.sc$delta=='down'])), annColors = list(osc=c('olivedrab1','purple4'), delta=c('navyblue', 'tomato3'),grh1=c('grey90', 'black')), annRow = list(grh1=Cao.scLS['WBGene00001707',]))


#(3)---matrisome---
#http://ec2-3-120-159-30.eu-central-1.compute.amazonaws.com:3838/ubuntu/ecm_analyzer/ 
#accessed on 26.10.2020
matrisome <- readRDS('Data/2020-10-26 Matrisome_results.rds')
matrisome <- data.frame(matrisome[,-c(1,5,8:16)])
rownames(matrisome) <- matrisome$WormBaseID
matrisome$class <- paste(matrisome$division, matrisome$category, sep='-')

genes.osc <- rownames(genesAnnot[genesAnnot$osc==T,])
genes.nonOsc <- rownames(genesAnnot[genesAnnot$osc==F,])
genes.up <- rownames(genesAnnot[genesAnnot$delta=='up',])
genes.down <- rownames(genesAnnot[genesAnnot$delta=='down',])

matrisome.osc <- matrisome[rownames(matrisome) %in% genes.osc,]
matrisome.nonOsc <- matrisome[rownames(matrisome) %in% genes.nonOsc,]
matrisome.up <- matrisome[rownames(matrisome) %in% genes.up,]
matrisome.down <- matrisome[rownames(matrisome) %in% genes.down,]

table(matrisome$class)
table(matrisome.osc$class)             
table(matrisome.nonOsc$class)
table(matrisome.up$class)
table(matrisome.down$class)
#-> enter manually into csv file

matrisome.comb <- data.frame(read.csv('Data/2020-10-26 Matrisome_results_all.csv', row.names = 1))
matrisome.comb.pred <- data.frame(pred.all=matrisome.comb$total/20000*length(genes))
matrisome.comb.pred$pred.osc <-matrisome.comb$total/20000*length(genes.osc)
matrisome.comb.pred$pred.nonOsc <-matrisome.comb$total/20000*length(genes.nonOsc)
matrisome.comb.pred$pred.up <-matrisome.comb$total/20000*length(genes.up)
matrisome.comb.pred$pred.down <-matrisome.comb$total/20000*length(genes.down)
rownames(matrisome.comb.pred) <- rownames(matrisome.comb)
matrisome.comb.fc <- (matrisome.comb[,-1]+8)/(matrisome.comb.pred+8)

par(mfrow=c(1,2))
plot.new()
barplot(t(matrisome.comb.fc-1),space=c(1,2),horiz = TRUE, las=TRUE,beside=TRUE,xaxt='n',col=brewer.pal(n = ncol(matrisome.comb.fc), name = "Greys"), xlab='Fold enrichment'); 
legend('topright', legend=colnames(matrisome.comb.fc), lty=1,col=brewer.pal(n = ncol(matrisome.comb.fc), name = "Greys"));
axis(1, at=seq(-1,2, by=0.5), labels=seq(-1,2, by=0.5)+1)  






