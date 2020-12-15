writeLines(capture.output(sessionInfo()), "grh-1_lowAuxinRange_sessionInfo.txt")

'
Analyze Luciferase Assay data Test64: grh-1::degron on low auxin concentration range

DESCRIPTION
-----------
This code allows the user to analyze the luciferase assay data of Test64 and generate figures of the GRH-1 manuscript

This script consists of the following steps: 
  (1) load data
  (2) calculate molt progression and plot
  (3) calculate stage durations
  (4) Quantify and plot M1-3 and I1-3 in animals that enter M4
  (5) quantification I1 and M1 in animals that do NOT show a M1 phenotype
  (6) determine the median luminescence levels in each molt of each worm


REFERENCES
----------
LuciferaseAssayAnalyzer (Meeuse, 2020)


USAGE
-----
Data format: 
file: .mat file from LuciferaseAssayAnalyzer containing molt and hatch annotations, luminescence values, and sample names
MoltProgr: .csv file containing mannual annotation of observed molts, i.e. observed molts were scored as 1, if molts were absent they were scored as 0


INPUT
-----
2018-04-09 low auxin range grh degron.mat
2019-02-28 low auxin range grh degron_Annotation_MoltPresence.csv

OUTPUT
------
survival plot: molt progression
boxplots: molt durations, luminescence intensity during molt
statistics: molt durations, luminescence intensity during molt


REQUIREMENTS
------------
This script requires the following packages to be installed:
R version 4.0.0 (2020-04-24)
RColorBrewer_1.1-2 
R.matlab_3.6.2 



"'


library('R.matlab')
library("RColorBrewer")

setwd("/tungstenfs/scratch/ggrossha/meeumilo/scripts/ggrosshans_grh-1_paper")


#(1)---Load data---
file <- readMat('Data/2018-04-09 low auxin range grh degron.mat')
valid <- as.factor(file$valid)
annotation <- unlist(file$samples)
annotation <- ordered(annotation, levels=c("ctr", "53 nM", "79 nM", "119 nM", "178 nM", "267 nM", '400 nM')) #order the annotation to get right order in boxplot
MoltProgr <- read.csv('Data/2019-02-28 low auxin range grh degron_Annotation_MoltPresence.csv', row.names=1)
MoltProgr$group <- annotation


#(2)---Molt progression and plotting---
MoltProgr_cum_ctr <- apply(MoltProgr[MoltProgr$group=='ctr',1:4],2,sum)
MoltProgr_frac_ctr <- (MoltProgr_cum_ctr/MoltProgr_cum_ctr[1])*100
MoltProgr_cum_53nM <- apply(MoltProgr[MoltProgr$group=='53 nM',1:4],2,sum)
MoltProgr_frac_53nM <- (MoltProgr_cum_53nM/MoltProgr_cum_53nM[1])*100
MoltProgr_cum_79nM <- apply(MoltProgr[MoltProgr$group=='79 nM',1:4],2,sum)
MoltProgr_frac_79nM <- (MoltProgr_cum_79nM/MoltProgr_cum_79nM[1])*100
MoltProgr_cum_119nM <- apply(MoltProgr[MoltProgr$group=='119 nM',1:4],2,sum)
MoltProgr_frac_119nM <- (MoltProgr_cum_119nM/MoltProgr_cum_119nM[1])*100
MoltProgr_cum_178nM <- apply(MoltProgr[MoltProgr$group=='178 nM',1:4],2,sum)
MoltProgr_frac_178nM <- (MoltProgr_cum_178nM/MoltProgr_cum_178nM[1])*100
MoltProgr_cum_267nM <- apply(MoltProgr[MoltProgr$group=='267 nM',1:4],2,sum)
MoltProgr_frac_267nM <- (MoltProgr_cum_267nM/MoltProgr_cum_267nM[1])*100
MoltProgr_cum_400nM <- apply(MoltProgr[MoltProgr$group=='400 nM',1:4],2,sum)
MoltProgr_frac_400nM <- (MoltProgr_cum_400nM/MoltProgr_cum_400nM[1])*100

MoltProgr_frac <- rbind(MoltProgr_frac_ctr, MoltProgr_frac_53nM, MoltProgr_frac_79nM, MoltProgr_frac_119nM, MoltProgr_frac_178nM, MoltProgr_frac_267nM, MoltProgr_frac_400nM)

col_step <- brewer.pal(n = 8, name = "Set1")

par(mfrow=c(3,4))
plot(stepfun(c(2,3,4), MoltProgr_frac[1,]), xlim=c(1,4), ylim=c(0,100), col=col_step[1], ylab='Molt completion (%)', xlab='Molt', main=paste('0 nM auxin ( n =', MoltProgr_cum_ctr[[1]],')')) 
plot(stepfun(c(2,3,4), MoltProgr_frac[2,]), xlim=c(1,4), ylim=c(0,100), col=col_step[2], ylab='Molt completion (%)', xlab='Molt', main=paste('53 nM auxin ( n =', MoltProgr_cum_53nM[[1]],')'))
plot(stepfun(c(2,3,4), MoltProgr_frac[3,]), xlim=c(1,4), ylim=c(0,100), col=col_step[3], ylab='Molt completion (%)', xlab='Molt', main=paste('79 nM auxin ( n =', MoltProgr_cum_79nM[[1]],')'))
plot(stepfun(c(2,3,4), MoltProgr_frac[4,]), xlim=c(1,4), ylim=c(0,100), col=col_step[4], ylab='Molt completion (%)', xlab='Molt', main=paste('119 nM auxin ( n =', MoltProgr_cum_119nM[[1]],')'))
plot(stepfun(c(2,3,4), MoltProgr_frac[5,]), xlim=c(1,4), ylim=c(0,100), col=col_step[5], ylab='Molt completion (%)', xlab='Molt', main=paste('178 nM auxin ( n =', MoltProgr_cum_178nM[[1]],')'))
plot(stepfun(c(2,3,4), MoltProgr_frac[6,]), xlim=c(1,4), ylim=c(0,100), col=col_step[7], ylab='Molt completion (%)', xlab='Molt', main=paste('267 nM auxin ( n =', MoltProgr_cum_267nM[[1]],')'))
plot(stepfun(c(2,3,4), MoltProgr_frac[7,]), xlim=c(1,4), ylim=c(0,100), col=col_step[8], ylab='Molt completion (%)', xlab='Molt', main=paste('400 nM auxin ( n =', MoltProgr_cum_400nM[[1]],')'))


#(3)---Stage durations---

M1 <- c()
M2 <- c()
M3 <- c()
M4 <- c()
I1 <- c()
I2 <- c()
I3 <- c()
I4 <- c()

for (i in 1:length(annotation)){
  M1[i] <- (file$Molt[1,2,i] - file$Molt[1,1,i])/6
  M2[i] <- (file$Molt[2,2,i] - file$Molt[2,1,i])/6
  M3[i] <- (file$Molt[3,2,i] - file$Molt[3,1,i])/6
  M4[i] <- (file$Molt[4,2,i] - file$Molt[4,1,i])/6
  I1[i] <- (file$Molt[1,1,i] - file$Hatch[1,i])/6
  I2[i] <- (file$Molt[2,1,i] - file$Molt[1,2,i])/6
  I3[i] <- (file$Molt[3,1,i] - file$Molt[2,2,i])/6
  I4[i] <- (file$Molt[4,1,i] - file$Molt[3,2,i])/6
}

df <- data.frame(condition = annotation, valid = valid, I1=I1,M1=M1,I2=I2,M2=M2,I3=I3,M3=M3,I4=I4, M4=M4)


#(4)---Quantify and plot M1-3 and I1-3 in animals that enter M4---
par(mfrow=c(2,3))
selection <- which(MoltProgr$M1==1 & MoltProgr$M2==1 &MoltProgr$M3==1 & MoltProgr$M4==1 &valid==1)
table(df$condition[selection])

boxplot(df$M1[selection] ~ df$condition[selection], data=df, ylab='M1 duration (h)', ylim=c(0.9*min(c(df$M1[selection],df$M2[selection],df$M3[selection])),1.1*max(c(df$M1[selection],df$M2[selection],df$M3[selection]))))
boxplot(df$M2[selection] ~ df$condition[selection], data=df, ylab='M2 duration (h)', ylim=c(0.9*min(c(df$M1[selection],df$M2[selection],df$M3[selection])),1.1*max(c(df$M1[selection],df$M2[selection],df$M3[selection]))))
boxplot(df$M3[selection] ~ df$condition[selection], data=df, ylab='M3 duration (h)', ylim=c(0.9*min(c(df$M1[selection],df$M2[selection],df$M3[selection])),1.1*max(c(df$M1[selection],df$M2[selection],df$M3[selection]))))

boxplot(df$I1[selection] ~ df$condition[selection], data=df, ylab='I1 duration (h)', ylim=c(0.9*min(c(df$I1[selection],df$I2[selection],df$I3[selection])),1.1*max(c(df$I1[selection],df$I2[selection],df$I3[selection]))))
boxplot(df$I2[selection] ~ df$condition[selection], data=df, ylab='I2 duration (h)', ylim=c(0.9*min(c(df$I1[selection],df$I2[selection],df$I3[selection])),1.1*max(c(df$I1[selection],df$I2[selection],df$I3[selection]))))
boxplot(df$I3[selection] ~ df$condition[selection], data=df, ylab='I3 duration (h)', ylim=c(0.9*min(c(df$I1[selection],df$I2[selection],df$I3[selection])),1.1*max(c(df$I1[selection],df$I2[selection],df$I3[selection]))))

#test normal distribution and equal variances
shapiro.test(df$M1[which(df$condition=='53 nM' & MoltProgr$M1==1 & MoltProgr$M2==1 &MoltProgr$M3==1 & MoltProgr$M4==1 &valid==1)]) #not normally distributed
var.test(df$M1[which(df$condition=='53 nM' & MoltProgr$M1==1 & MoltProgr$M2==1 &MoltProgr$M3==1 & MoltProgr$M4==1 &valid==1)],df$M1[which(df$condition=='ctr' & MoltProgr$M1==1 & MoltProgr$M2==1 &MoltProgr$M3==1 & MoltProgr$M4==1 &valid==1)])

selection_ctl <- selection[selection %in% which(df$condition=='ctr')]
selection_53nM <- selection[selection %in% which(df$condition=='53 nM')]
selection_79nM <- selection[selection %in% which(df$condition=='79 nM')]
selection_119nM <- selection[selection %in% which(df$condition=='119 nM')]

pval_M1_53nM <- wilcox.test(df$M1[selection_53nM], df$M1[selection_ctl])$p.value
pval_M1_79nM <- wilcox.test(df$M1[selection_79nM], df$M1[selection_ctl])$p.value
pval_M1_119nM <- wilcox.test(df$M1[selection_119nM], df$M1[selection_ctl])$p.value
pval_M2_53nM <- wilcox.test(df$M2[selection_53nM], df$M2[selection_ctl])$p.value
pval_M2_79nM <- wilcox.test(df$M2[selection_79nM], df$M2[selection_ctl])$p.value
pval_M2_119nM <- wilcox.test(df$M2[selection_119nM], df$M2[selection_ctl])$p.value
pval_M3_53nM <- wilcox.test(df$M3[selection_53nM], df$M3[selection_ctl])$p.value
pval_M3_79nM <- wilcox.test(df$M3[selection_79nM], df$M3[selection_ctl])$p.value
pval_M3_119nM <- wilcox.test(df$M3[selection_119nM], df$M3[selection_ctl])$p.value
pval_M <- rbind(pval_M1_53nM,pval_M1_79nM,pval_M1_119nM,pval_M2_53nM,pval_M2_79nM,pval_M2_119nM,pval_M3_53nM,pval_M3_79nM,pval_M3_119nM)

pval_I1_53nM <- wilcox.test(df$I1[selection_53nM], df$I1[selection_ctl])$p.value
pval_I1_79nM <- wilcox.test(df$I1[selection_79nM], df$I1[selection_ctl])$p.value
pval_I1_119nM <- wilcox.test(df$I1[selection_119nM], df$I1[selection_ctl])$p.value
pval_I2_53nM <- wilcox.test(df$I2[selection_53nM], df$I2[selection_ctl])$p.value
pval_I2_79nM <- wilcox.test(df$I2[selection_79nM], df$I2[selection_ctl])$p.value
pval_I2_119nM <- wilcox.test(df$I2[selection_119nM], df$I2[selection_ctl])$p.value
pval_I3_53nM <- wilcox.test(df$I3[selection_53nM], df$I3[selection_ctl])$p.value
pval_I3_79nM <- wilcox.test(df$I3[selection_79nM], df$I3[selection_ctl])$p.value
pval_I3_119nM <- wilcox.test(df$I3[selection_119nM], df$I3[selection_ctl])$p.value
pval_I <- rbind(pval_I1_53nM,pval_I1_79nM,pval_I1_119nM,pval_I2_53nM,pval_I2_79nM,pval_I2_119nM,pval_I3_53nM,pval_I3_79nM,pval_I3_119nM)

pval <- cbind(pval_M, pval_I)
colnames(pval) <- c('M1', 'I1')
write.csv(pval, 'Plots/2020-06-17 pval_Wilcox_M1-3_I1-3.csv')


#(5)---quantification I1 and M1 in animals that do NOT show a M1 phenotype---
selection2 <- which(MoltProgr$M1==1 & MoltProgr$M2==1 &valid==1)
table(df$condition[selection2])
par(mfrow=c(1,1))
boxplot(df$M1[selection2] ~ df$condition[selection2], data=df, ylab='M1 duration (h)')
boxplot(df$I1[selection2] ~ df$condition[selection2], data=df, ylab='I1 duration (h)')

pval_M1_53nM_normal <- wilcox.test(df$M1[which(df$condition=='53 nM' & MoltProgr$M1==1 & MoltProgr$M2==1 &valid==1)], df$M1[which(df$condition=='ctr' & MoltProgr$M1==1 & MoltProgr$M2==1 &valid==1)])$p.value
pval_M1_79nM_normal <- wilcox.test(df$M1[which(df$condition=='79 nM' & MoltProgr$M1==1 & MoltProgr$M2==1 &valid==1)], df$M1[which(df$condition=='ctr' & MoltProgr$M1==1 & MoltProgr$M2==1 &valid==1)])$p.value
pval_M1_119nM_normal <- wilcox.test(df$M1[which(df$condition=='119 nM' & MoltProgr$M1==1 & MoltProgr$M2==1 &valid==1)], df$M1[which(df$condition=='ctr' & MoltProgr$M1==1 & MoltProgr$M2==1 &valid==1)])$p.value
pval_M1_178nM_normal <- wilcox.test(df$M1[which(df$condition=='178 nM' & MoltProgr$M1==1 & MoltProgr$M2==1 &valid==1)], df$M1[which(df$condition=='ctr' & MoltProgr$M1==1 & MoltProgr$M2==1 &valid==1)])$p.value
pval_M1_267nM_normal <- wilcox.test(df$M1[which(df$condition=='267 nM' & MoltProgr$M1==1 & MoltProgr$M2==1 &valid==1)], df$M1[which(df$condition=='ctr' & MoltProgr$M1==1 & MoltProgr$M2==1 &valid==1)])$p.value
pval_M1_normal <- rbind(pval_M1_53nM_normal, pval_M1_79nM_normal, pval_M1_119nM_normal, pval_M1_178nM_normal, pval_M1_267nM_normal)

pval_I1_53nM_normal <- wilcox.test(df$I1[which(df$condition=='53 nM' & MoltProgr$M1==1 & MoltProgr$M2==1 &valid==1)], df$I1[which(df$condition=='ctr' & MoltProgr$M1==1 & MoltProgr$M2==1 &valid==1)])$p.value
pval_I1_79nM_normal <- wilcox.test(df$I1[which(df$condition=='79 nM' & MoltProgr$M1==1 & MoltProgr$M2==1 &valid==1)], df$I1[which(df$condition=='ctr' & MoltProgr$M1==1 & MoltProgr$M2==1 &valid==1)])$p.value
pval_I1_119nM_normal <- wilcox.test(df$I1[which(df$condition=='119 nM' & MoltProgr$M1==1 & MoltProgr$M2==1 &valid==1)], df$I1[which(df$condition=='ctr' & MoltProgr$M1==1 & MoltProgr$M2==1 &valid==1)])$p.value
pval_I1_178nM_normal <- wilcox.test(df$I1[which(df$condition=='178 nM' & MoltProgr$M1==1 & MoltProgr$M2==1 &valid==1)], df$I1[which(df$condition=='ctr' & MoltProgr$M1==1 & MoltProgr$M2==1 &valid==1)])$p.value
pval_I1_267nM_normal<- wilcox.test(df$I1[which(df$condition=='267 nM' & MoltProgr$M1==1 & MoltProgr$M2==1 &valid==1)], df$I1[which(df$condition=='ctr' & MoltProgr$M1==1 & MoltProgr$M2==1 &valid==1)])$p.value
pval_I1_normal <- rbind(pval_I1_53nM_normal, pval_I1_79nM_normal, pval_I1_119nM_normal, pval_I1_178nM_normal, pval_I1_267nM_normal)

pval_L1_normal <- cbind(pval_M1_normal, pval_I1_normal)
colnames(pval_L1_normal) <- c('M1', 'I1')
write.csv(pval_L1_normal, 'Plots/2020-06-17 pval_L1_normal.csv')

# quantification of I1 in 400 nM
boxplot(df$I1[which(df$condition=='ctr' & valid==1)], df$I1[which(df$condition=='267 nM' & valid==1 &MoltProgr$M2==0)], df$I1[which(df$condition=='400 nM' & valid==1)], xlab=c('veh', '267 nM','400 nM'), ylab='I1 duration (h)') 
table(c(df$condition[which(df$condition=='ctr' & valid==1)], df$condition[which(df$condition=='267 nM' & valid==1 &MoltProgr$M2==0)], df$condition[which(df$condition=='400 nM' & valid==1)]))
wilcox.test(df$I1[which(df$condition=='ctr' & valid==1)], df$I1[which(df$condition=='400 nM' & valid==1)]) 
wilcox.test(df$I1[which(df$condition=='ctr' & valid==1)], df$I1[which(df$condition=='267 nM' & valid==1 &MoltProgr$M2==0)])


#(6)---determine the median luminescence levels in each molt of each worm---
Lumi <- t(file$X)
M1_lumi <- c()
M2_lumi <- c()
M3_lumi <- c()
M4_lumi <- c()

for (i in 1:nrow(Lumi)){
  M1_lumi[i] <- median(Lumi[i,file$Molt[1,1,i]:file$Molt[1,2,i]]) #median because less affected by wrongly set points
  M2_lumi[i] <- median(Lumi[i,file$Molt[2,1,i]:file$Molt[2,2,i]])
  M3_lumi[i] <- median(Lumi[i,file$Molt[3,1,i]:file$Molt[3,2,i]])
  M4_lumi[i] <- median(Lumi[i,file$Molt[4,1,i]:file$Molt[4,2,i]])
}

df_lumi <- df[,1:2]
df_lumi$M1_lumi <- M1_lumi
df_lumi$M2_lumi <- M2_lumi
df_lumi$M3_lumi <- M3_lumi
df_lumi$M4_lumi <- M4_lumi

df_lumi_sel <- df_lumi[selection,]

M1_lumi_ctr_median <- median(df_lumi_sel$M1_lumi[which(df_lumi_sel$condition=='ctr')])
M2_lumi_ctr_median <- median(df_lumi_sel$M2_lumi[which(df_lumi_sel$condition=='ctr')])
M3_lumi_ctr_median <- median(df_lumi_sel$M3_lumi[which(df_lumi_sel$condition=='ctr')])
M4_lumi_ctr_median <- median(df_lumi_sel$M4_lumi[which(df_lumi_sel$condition=='ctr')])

par(mfrow=c(2,4))
boxplot(df_lumi_sel$M1_lumi/M1_lumi_ctr_median ~ df_lumi_sel$condition, data=df_lumi_sel, main='Median Luminescence Molt 1', ylab='fold change luminescence', ylim=c(0,5))
boxplot(df_lumi_sel$M2_lumi/M2_lumi_ctr_median ~ df_lumi_sel$condition, data=df_lumi_sel, main='Median Luminescence Molt 2', ylab='fold change luminescence', ylim=c(0,5))
boxplot(df_lumi_sel$M3_lumi/M3_lumi_ctr_median ~ df_lumi_sel$condition, data=df_lumi_sel, main='Median Luminescence Molt 3', ylab='fold change luminescence', ylim=c(0,15))
boxplot(df_lumi_sel$M4_lumi/M4_lumi_ctr_median ~ df_lumi_sel$condition, data=df_lumi_sel, main='Median Luminescence Molt 4', ylab='fold change luminescence', ylim=c(0,15))

shapiro.test(df_lumi_sel$M4_lumi[which(df_lumi_sel$condition=='53 nM')])
var.test(df_lumi_sel$M4_lumi[which(df_lumi_sel$condition=='53 nM')], df_lumi_sel$M4_lumi[which(df_lumi_sel$condition=='ctr')])$p.value
wilcox.test(df_lumi_sel$M4_lumi[which(df_lumi_sel$condition=='53 nM')], df_lumi_sel$M4_lumi[which(df_lumi_sel$condition=='ctr')])$p.value

df_lumi_sel_53nM <- df_lumi_sel[which(df_lumi_sel$condition=='53 nM' | df_lumi_sel$condition=='ctr'),]
df_lumi_sel_79nM <- df_lumi_sel[which(df_lumi_sel$condition=='79 nM' | df_lumi_sel$condition=='ctr'),]
df_lumi_sel_119nM <- df_lumi_sel[which(df_lumi_sel$condition=='119 nM' | df_lumi_sel$condition=='ctr'),]
pval_lumi_53nM <- sapply(df_lumi_sel_53nM[,3:6], function(i) {wilcox.test(i ~ df_lumi_sel_53nM$condition)$p.value})
pval_lumi_79nM <- sapply(df_lumi_sel_79nM[,3:6], function(i) {wilcox.test(i ~ df_lumi_sel_79nM$condition)$p.value})
pval_lumi_119nM <- sapply(df_lumi_sel_119nM[,3:6], function(i) {wilcox.test(i ~ df_lumi_sel_119nM$condition)$p.value})
pval_lumi <- data.frame(p_53nM = pval_lumi_53nM, p_79nM = pval_lumi_79nM, p_119nM = pval_lumi_119nM)
write.csv(pval_lumi, '2020-06-17 pval_lumi_53nM_79nM_119nM.csv')









