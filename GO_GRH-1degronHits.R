writeLines(capture.output(date(),sessionInfo()), "GO_GRH-1degronHits_sessionInfo.txt")
library('ggplot2')

'
GO terms to barplot 

DESCRIPTION
-----------
This code allows the user to convert the FoldEnrichment and the FalseDiscoveryRate of the GO terms from PANTHER to a barplot
This script consists of the following steps: 
(1) load data
(2) set the order of the GO terms 
(3) make the barplot


REFERENCES
----------
GO term enrichment calculated using: http://www.pantherdb.org/

    
INPUT
-----
Data format: csv file with columns: GO, FoldEnr, Pval, FDR
settings: accessed on 27.10.2020, PANTHER Overrepresentation Test released 20200728
          analyzed list: 2020-10-26 DiffExpr_RNAseq1779.csv (623 genes mapped, 10 unmapped)
          reference list: C.elegans genome (19961 genes)
          test type: Fisher exact with FDR correction
          

OUTPUT
------
Results are stored as: barplot 
The output file has the following information: with every row one GO term, the x-axis=log2(FoldEnr), and the color=log10(FDR)


REQUIREMENTS
------------
ggplot2_3.3.2


"'
setwd("/tungstenfs/scratch/ggrossha/meeumilo/scripts/ggrosshans_grh-1_paper")
library(ggplot2)

#---load data---

GO_BP <- read.csv('Data/2020-10-27 GOterm_BiolProc.csv')
GO_BP <- GO_BP[GO_BP$Fold.Enrichment>1,]
GO_MF <- read.csv('Data/2020-10-27 GOterm_MolFun.csv')
GO_MF <- GO_MF[GO_MF$Fold.Enrichment>1,]

  # reverse the order of GO such that the order of the barplot labels goes from high to low FoldEnr
  GO_BP$GO <- factor(GO_BP$GO, levels = rev(GO_BP$GO)) 
  GO_MF$GO <- factor(GO_MF$GO, levels = rev(GO_MF$GO)) 

#---make barplots---
  
  # Determine the ranges of the color scheme and the x-axis
  lim_y <- c(0,max(log2(range(c(GO_BP$Fold.Enrichment,GO_MF$Fold.Enrichment))))) #y-axis will become x-axis by coord_flip()

  # Make barplot with individual color scheme
  p1<-ggplot(data=GO_BP, aes(x=GO, y=log2(Fold.Enrichment), fill=log10(FDR))) +
    geom_bar(stat="identity") + labs(title='GO BioProc') + scale_fill_gradient2(low="black", mid="grey", high="white")
  p1 + coord_flip() + scale_y_continuous(limits=lim_y) 
  
  p2<-ggplot(data=GO_MF, aes(x=GO, y=log2(Fold.Enrichment), fill=log10(FDR))) +
    geom_bar(stat="identity") + labs(title='GO MolFun')  + scale_fill_gradient2(low="black", mid="grey", high="white")  
  p2 + coord_flip() + scale_y_continuous(limits=lim_y) 



