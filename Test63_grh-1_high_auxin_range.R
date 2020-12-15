writeLines(capture.output(sessionInfo()), "Test63_grh-1_high_auxin_range_sessionInfo.txt")

'
Analyze Luciferase Assay data Test63: grh-1::degron on high auxin concentration range

DESCRIPTION
-----------
This code allows the user to analyze the luciferase assay data of Test64 and generate figures of the GRH-1 manuscript

This script consists of the following steps: 
  (1) load data
  (2) calculate molt progression and plot
  (3) calculate stage durations
  (4) Quantify and plot M1 and I1 in animals that do NOT show M1 phenotype
  (5) Quantify and plot I1 in animals that complete M1


REFERENCES
----------
LuciferaseAssayAnalyzer (Meeuse et al, MSB, 2020)


USAGE
-----
Data format: 
file: .mat file from LuciferaseAssayAnalyzer containing molt and hatch annotations, luminescence values, and sample names
MoltProgr: .csv file containing mannual annotation of observed molts, i.e. observed molts were scored as 1, if molts were absent they were scored as 0


INPUT
-----
2018-03-27 grh degron auxin range.mat
2020-06-24 Annotation_MoltPresence.csv

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


setwd("/tungstenfs/scratch/ggrossha/meeumilo/scripts/ggrosshans_grh-1_paper")
library('R.matlab')
library("RColorBrewer")



#(1)---Load data---
file <- readMat('Data/2018-03-27 grh-1_highAuxRange.mat')
valid <- as.factor(file$valid)
annotation <- unlist(file$samples)
annotation <- ordered(annotation, levels=c("ctr", "60.9 nM", "243 nM", "975 nM", "3.9 uM", "15.6 uM", '62.5 uM', '250 uM', '1 mM')) #order the annotation to get right order in boxplot
MoltProgr <- read.csv('Data/2020-06-24 grh-1_highAuxRange_Annotation_MoltPresence.csv', row.names=1)
MoltProgr$group <- annotation
MoltProgr <- MoltProgr[valid==1,] #some animals hatch (valid==1) but do not molt (M1==0)

#(2)---Molt progression and plotting---
MoltProgr_cum_ctr <- apply(MoltProgr[MoltProgr$group=='ctr',1:4],2,sum)
MoltProgr_frac_ctr <- (MoltProgr_cum_ctr/nrow(MoltProgr[MoltProgr$group=='ctr',]))*100
MoltProgr_cum_60nM <- apply(MoltProgr[MoltProgr$group=='60.9 nM',1:4],2,sum)
MoltProgr_frac_60nM <- (MoltProgr_cum_60nM/nrow(MoltProgr[MoltProgr$group=='60.9 nM',]))*100
MoltProgr_cum_243nM <- apply(MoltProgr[MoltProgr$group=='243 nM',1:4],2,sum)
MoltProgr_frac_243nM <- (MoltProgr_cum_243nM/nrow(MoltProgr[MoltProgr$group=='243 nM',]))*100
MoltProgr_cum_975nM <- apply(MoltProgr[MoltProgr$group=='975 nM',1:4],2,sum)
MoltProgr_frac_975nM <- (MoltProgr_cum_975nM/nrow(MoltProgr[MoltProgr$group=='975 nM',]))*100
MoltProgr_cum_4uM <- apply(MoltProgr[MoltProgr$group=='3.9 uM',1:4],2,sum)
MoltProgr_frac_4uM <- (MoltProgr_cum_4uM/nrow(MoltProgr[MoltProgr$group=='3.9 uM',]))*100
MoltProgr_cum_15uM <- apply(MoltProgr[MoltProgr$group=='15.6 uM',1:4],2,sum)
MoltProgr_frac_15uM <- (MoltProgr_cum_15uM/nrow(MoltProgr[MoltProgr$group=='15.6 uM',]))*100
MoltProgr_cum_62uM <- apply(MoltProgr[MoltProgr$group=='62.5 uM',1:4],2,sum)
MoltProgr_frac_62uM <- (MoltProgr_cum_62uM/nrow(MoltProgr[MoltProgr$group=='62.5 uM',]))*100
MoltProgr_cum_250uM <- apply(MoltProgr[MoltProgr$group=='250 uM',1:4],2,sum)
MoltProgr_frac_250uM <- (MoltProgr_cum_250uM/nrow(MoltProgr[MoltProgr$group=='250 uM',]))*100
MoltProgr_cum_1mM <- apply(MoltProgr[MoltProgr$group=='1 mM',1:4],2,sum)
MoltProgr_frac_1mM <- (MoltProgr_cum_1mM/nrow(MoltProgr[MoltProgr$group=='1 mM',]))*100

MoltProgr_frac <- rbind(MoltProgr_frac_ctr, MoltProgr_frac_60nM, MoltProgr_frac_243nM, MoltProgr_frac_975nM, MoltProgr_frac_4uM, MoltProgr_frac_15uM, MoltProgr_frac_62uM,MoltProgr_frac_250uM,MoltProgr_frac_1mM)

col_step <- brewer.pal(n = 9, name = "Set3")

par(mfrow=c(3,4))
plot(stepfun(c(2,3,4), MoltProgr_frac[1,]), xlim=c(1,4), ylim=c(0,100), col=col_step[1], ylab='Molt completion (%)', xlab='Molt', main=paste('0 nM auxin ( n =', nrow(MoltProgr[MoltProgr$group=='ctr',]),')')) 
plot(stepfun(c(2,3,4), MoltProgr_frac[2,]), xlim=c(1,4), ylim=c(0,100), col=col_step[2], ylab='Molt completion (%)', xlab='Molt', main=paste('60.9 nM auxin ( n =', nrow(MoltProgr[MoltProgr$group=='60.9 nM',]),')'))
plot(stepfun(c(2,3,4), MoltProgr_frac[3,]), xlim=c(1,4), ylim=c(0,100), col=col_step[3], ylab='Molt completion (%)', xlab='Molt', main=paste('243 nM auxin ( n =', nrow(MoltProgr[MoltProgr$group=='243 nM',]),')'))
plot(stepfun(c(2,3,4), MoltProgr_frac[4,]), xlim=c(1,4), ylim=c(0,100), col=col_step[4], ylab='Molt completion (%)', xlab='Molt', main=paste('975 nM auxin ( n =', nrow(MoltProgr[MoltProgr$group=='975 nM',]),')'))
plot(stepfun(c(2,3,4), MoltProgr_frac[5,]), xlim=c(1,4), ylim=c(0,100), col=col_step[5], ylab='Molt completion (%)', xlab='Molt', main=paste('3.9 uM auxin ( n =', nrow(MoltProgr[MoltProgr$group=='3.9 uM',]),')'))
plot(stepfun(c(2,3,4), MoltProgr_frac[6,]), xlim=c(1,4), ylim=c(0,100), col=col_step[6], ylab='Molt completion (%)', xlab='Molt', main=paste('15.6 uM auxin ( n =', nrow(MoltProgr[MoltProgr$group=='15.6 uM',]),')'))
plot(stepfun(c(2,3,4), MoltProgr_frac[7,]), xlim=c(1,4), ylim=c(0,100), col=col_step[7], ylab='Molt completion (%)', xlab='Molt', main=paste('62.5 uM auxin ( n =', nrow(MoltProgr[MoltProgr$group=='62.5 uM',]),')'))
plot(stepfun(c(2,3,4), MoltProgr_frac[8,]), xlim=c(1,4), ylim=c(0,100), col=col_step[8], ylab='Molt completion (%)', xlab='Molt', main=paste('250 uM auxin ( n =', nrow(MoltProgr[MoltProgr$group=='250 uM',]),')'))
plot(stepfun(c(2,3,4), MoltProgr_frac[9,]), xlim=c(1,4), ylim=c(0,100), col=col_step[9], ylab='Molt completion (%)', xlab='Molt', main=paste('1 mM auxin ( n =', nrow(MoltProgr[MoltProgr$group=='1 mM',]),')'))



#(3)---Stage durations---
M1 <- c()
M2 <- c()
M3 <- c()
M4 <- c()
I1 <- c()
I2 <- c()
I3 <- c()
I4 <- c()
M1_in <- c()

for (i in 1:length(annotation)){
  M1[i] <- (file$Molt[1,2,i] - file$Molt[1,1,i])/6
  M2[i] <- (file$Molt[2,2,i] - file$Molt[2,1,i])/6
  M3[i] <- (file$Molt[3,2,i] - file$Molt[3,1,i])/6
  M4[i] <- (file$Molt[4,2,i] - file$Molt[4,1,i])/6
  I1[i] <- (file$Molt[1,1,i] - file$Hatch[1,i])/6
  I2[i] <- (file$Molt[2,1,i] - file$Molt[1,2,i])/6
  I3[i] <- (file$Molt[3,1,i] - file$Molt[2,2,i])/6
  I4[i] <- (file$Molt[4,1,i] - file$Molt[3,2,i])/6
  M1_in[i] <- file$Molt[1,1,i]/6
}

df <- data.frame(condition = annotation, valid = valid, M1_in=M1_in, I1=I1,M1=M1,I2=I2,M2=M2,I3=I3,M3=M3,I4=I4, M4=M4)
df <- df[valid==1,]


#(4)---Quantify and plot M1 and I1 in animals that do NOT show M1 phenotype---
df_M1 <- df[which(MoltProgr$M1==1 & MoltProgr$M2==1),] #NOT show M1 phenotype = worms that molt into M2
par(mfrow=c(2,3))
table(df_M1$condition)
boxplot(M1 ~ condition, data=df_M1, las=2, main = '', ylab = 'Molt 1 duration (h)', xlab='')
boxplot(I1 ~ condition, data=df_M1, las=2, main = '', ylab = 'Intermolt 1 duration (h)', xlab='')

#-statistics
shapiro.test(df_M1$M1[which(df_M1$condition=='ctr')]) #not normally distributed
var.test(df_M1$M1[which(df_M1$condition=='ctr')],df_M1$M1[which(df_M1$condition=='60.9 nM')])

pval_M1_60nM <- wilcox.test(df_M1$M1[which(df_M1$condition=='60.9 nM')], df_M1$M1[which(df_M1$condition=='ctr')])$p.value
pval_M1_243nM <- wilcox.test(df_M1$M1[which(df_M1$condition=='243 nM')], df_M1$M1[which(df_M1$condition=='ctr')])$p.value
pval_M1 <- rbind(pval_M1_60nM, pval_M1_243nM)
pval_I1_60nM <- wilcox.test(df_M1$I1[which(df_M1$condition=='60.9 nM')], df_M1$I1[which(df_M1$condition=='ctr')])$p.value
pval_I1_243nM <- wilcox.test(df_M1$I1[which(df_M1$condition=='243 nM')], df_M1$I1[which(df_M1$condition=='ctr')])$p.value
pval_I1 <- rbind(pval_I1_60nM, pval_I1_243nM)

pval_L1 <- cbind(pval_M1, pval_I1)
colnames(pval_L1) <- c('M1', 'I1')
write.csv(pval_L1, 'Plots/LucAssays/Test63_grh-1_highAuxRange/2020-06-24 pval_L1.csv')


#(5)---Quantify and plot I1 in all animals---
df_M1all <- df[which(MoltProgr$M1==1),] #all animals that complete M1, they may show phenotype, but also may not
table(df_M1all$condition)
boxplot(I1 ~ condition, data=df_M1all, las=2, main = '', ylab = 'Intermolt 1 duration (h)', xlab='')

#-statistics
pval_I1all_60nM <- wilcox.test(df_M1all$I1[which(df_M1all$condition=='60.9 nM')], df_M1all$I1[which(df_M1all$condition=='ctr')])$p.value
pval_I1all_243nM <- wilcox.test(df_M1all$I1[which(df_M1all$condition=='243 nM')], df_M1all$I1[which(df_M1all$condition=='ctr')])$p.value
pval_I1all_975nM <- wilcox.test(df_M1all$I1[which(df_M1all$condition=='975 nM')], df_M1all$I1[which(df_M1all$condition=='ctr')])$p.value
pval_I1all_4uM <- wilcox.test(df_M1all$I1[which(df_M1all$condition=='3.9 uM')], df_M1all$I1[which(df_M1all$condition=='ctr')])$p.value
pval_I1all_15uM <- wilcox.test(df_M1all$I1[which(df_M1all$condition=='15.6 uM')], df_M1all$I1[which(df_M1all$condition=='ctr')])$p.value
pval_I1all_62uM <- wilcox.test(df_M1all$I1[which(df_M1all$condition=='62.5 uM')], df_M1all$I1[which(df_M1all$condition=='ctr')])$p.value
pval_I1all_250uM <- wilcox.test(df_M1all$I1[which(df_M1all$condition=='250 uM')], df_M1all$I1[which(df_M1all$condition=='ctr')])$p.value
pval_I1all_1mM <- wilcox.test(df_M1all$I1[which(df_M1all$condition=='1 mM')], df_M1all$I1[which(df_M1all$condition=='ctr')])$p.value


pval_I1all <- rbind(pval_I1all_60nM, pval_I1all_243nM, pval_I1all_975nM, pval_I1all_4uM, pval_I1all_15uM, pval_I1all_62uM, pval_I1all_250uM, pval_I1all_1mM)
write.csv(pval_I1all, 'Plots/LucAssays/Test63_grh-1_highAuxRange/2020-06-24 pval_I1_all.csv')








