writeLines(capture.output(date(),sessionInfo()), 'Hatching_GRH-1mutants_sessioninfo.txt')

'
Hatching experiment

DESCRIPTION
-----------
Analysis was performed by: Milou Meeuse

This code allows the user to make barplots 
  

EXPERIMENT DESCRIPTION
----------------------
Experiment was performed by: Anca Neagu
2020-11-20: • pick 10 L4s • let them develop to egg-laying adults for ~2days • pick 60 eggs on two 3.5cm plates and let them hatch for >24h
2020-11-27: Let gravid adults lay eggs for 3h. Pick 100 eggs on fresh plates. Check hatching after 24h. (in two batches coming from the same animals)

animals:  grh-1::gfp (named as grh)
          grh-1(dDseg2)::gfp/hT2 (named as grh_Dseg2)
          grh-1(0)::hT2 (named as grh_null)
          

REFERENCES
----------


INPUT
-----
- table with % of animals 


OUTPUT
------
- barplot


REQUIREMENTS
------------
 ggplot2_3.3.2


"'

setwd("/tungstenfs/scratch/ggrossha/meeumilo/scripts/ggrosshans_grh-1_paper")
library(ggplot2)

hatch <- read.csv('Data/2020-12-08 HatchingExp_Table.csv')
hatch.mean <- aggregate(hatch$percentage,list(hatch$genotype,hatch$phenotype),mean)
colnames(hatch.mean) <- colnames(hatch)
hatch.mean$sd <- aggregate(hatch$percentage,list(hatch$genotype,hatch$phenotype),sd)$x

ggplot(data=hatch.mean, aes(x=genotype, y=percentage, fill=phenotype)) + 
  geom_bar(position = "dodge", stat = "identity") + theme_classic()+
  geom_errorbar( aes( ymin=percentage-sd, ymax=percentage+sd), position=position_dodge(width=0.9),width=0.5, colour="orange", size=0.5)+
  scale_fill_manual("legend", values = c("grey90", "grey50","black"))


