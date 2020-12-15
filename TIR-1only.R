writeLines(capture.output(sessionInfo()), "TIR-1only_sessionInfo.txt")

library('R.matlab')
library('reshape2')
library('ggplot2')
library('survminer')
library('stats')
library("RColorBrewer")
setwd("/tungstenfs/scratch/ggrossha/meeumilo/scripts/ggrosshans_grh-1_paper")

file <- readMat('Data/2016-11-01 Test42_TIR1only.mat') #worms with 4 molts are annotated
valid <- as.factor(file$valid)
annotation <- unlist(file$samples)
annotation <- ordered(annotation, levels=c("pMM002.4 - EV RNAi", "pMM002.4 - nhr-25 RNAi",'pMM002.4 x CA1200 - 0.25% eth','pMM002.4 x CA1200 - 1% eth',"pMM002.4 x CA1200 - 3.9 uM auxin ",'pMM002.4 x CA1200 - 15.6 uM auxin ','pMM002.4 x CA1200 - 62.5 uM auxin ','pMM002.4 x CA1200 - 250 uM auxin ','pMM002.4 x CA1200 - 1 mM auxin ','pMM002.4 x CA1200 - 4 mM auxin ','pMM002.4 x KRY85 - 0.25% eth','pMM002.4 x KRY85 - 1% eth',"pMM002.4 x KRY85 - 3.9 uM auxin ",'pMM002.4 x KRY85 - 15.6 uM auxin ','pMM002.4 x KRY85 - 62.5 uM auxin ','pMM002.4 x KRY85 - 250 uM auxin ','pMM002.4 x KRY85 - 1 mM auxin ','pMM002.4 x KRY85 - 4 mM auxin '))
moltCompletion <- read.csv('Data/MoltCompletion_Test42_TIR1only.csv', row.names = 1) #E is excluded mostly because of empty well


##- duration stages
M1 <- c()
M2 <- c()
M3 <- c()
M4 <- c()
I1 <- c()
I2 <- c()
I3 <- c()
I4 <- c()
L1 <- c()
L2 <- c()
L3 <- c()
L4 <- c()

for (i in 1:length(valid)){
  M1[i] <- (file$Molt[1,2,i] - file$Molt[1,1,i])/6
  M2[i] <- (file$Molt[2,2,i] - file$Molt[2,1,i])/6
  M3[i] <- (file$Molt[3,2,i] - file$Molt[3,1,i])/6
  M4[i] <- (file$Molt[4,2,i] - file$Molt[4,1,i])/6
  I1[i] <- (file$Molt[1,1,i] - file$Hatch[1,i])/6
  I2[i] <- (file$Molt[2,1,i] - file$Molt[1,2,i])/6
  I3[i] <- (file$Molt[3,1,i] - file$Molt[2,2,i])/6
  I4[i] <- (file$Molt[4,1,i] - file$Molt[3,2,i])/6
  L1[i] <- (file$Molt[1,2,i] - file$Hatch[1,i])/6
  L2[i] <- (file$Molt[2,2,i] - file$Molt[1,2,i])/6
  L3[i] <- (file$Molt[3,2,i] - file$Molt[2,2,i])/6
  L4[i] <- (file$Molt[4,2,i] - file$Molt[3,2,i])/6
}

df <- data.frame(condition = annotation, valid = valid, M1=M1, M2=M2, M3=M3, M4=M4, I1=I1, I2=I2, I3=I3, I4=I4, L1=L1, L2=L2, L3=L3, L4=L4)
df <- df[df$valid==1,]

#quantification Tir-1 only: effect of auxin on development
df_TIR <- df[(grepl('CA1200', df$condition)),]
df_TIR$condition <- droplevels(df_TIR$condition, exclude=unique(df$condition[!(grepl('CA1200', df$condition))]))

par(mfrow=c(2,4))
boxplot(df_TIR$M1 ~ df_TIR$condition, data=df_TIR, main='Molt 1', las=2)
boxplot(df_TIR$M2 ~ df_TIR$condition, data=df_TIR, main='Molt 2', las=2)
boxplot(df_TIR$M3 ~ df_TIR$condition, data=df_TIR, main='Molt 3', las=2)
boxplot(df_TIR$M4 ~ df_TIR$condition, data=df_TIR, main='Molt 4', las=2)
par(mfrow=c(2,4))
boxplot(df_TIR$I1 ~ df_TIR$condition, data=df_TIR, main='Intermolt 1', las=2)
boxplot(df_TIR$I2 ~ df_TIR$condition, data=df_TIR, main='Intermolt 2', las=2)
boxplot(df_TIR$I3 ~ df_TIR$condition, data=df_TIR, main='Intermolt 3', las=2)
boxplot(df_TIR$I4 ~ df_TIR$condition, data=df_TIR, main='Intermolt 4', las=2)

wilcox.test(df_TIR$M4[which(df_TIR$condition=='pMM002.4 x CA1200 - 1% eth')], df_TIR$M4[which(df_TIR$condition=='pMM002.4 x CA1200 - 4 mM auxin ')])$p.value
pval_4mM <- sapply(df_TIR[which(df_TIR$condition=='pMM002.4 x CA1200 - 1% eth' | df_TIR$condition=='pMM002.4 x CA1200 - 4 mM auxin '),3:10], function(i) {wilcox.test(i ~ df_TIR$condition[which(df_TIR$condition=='pMM002.4 x CA1200 - 1% eth' | df_TIR$condition=='pMM002.4 x CA1200 - 4 mM auxin ')])$p.value})
pval_1mM <- sapply(df_TIR[which(df_TIR$condition=='pMM002.4 x CA1200 - 0.25% eth' | df_TIR$condition=='pMM002.4 x CA1200 - 1 mM auxin '),3:10], function(i) {wilcox.test(i ~ df_TIR$condition[which(df_TIR$condition=='pMM002.4 x CA1200 - 0.25% eth' | df_TIR$condition=='pMM002.4 x CA1200 - 1 mM auxin ')])$p.value})
pval_250uM <- sapply(df_TIR[which(df_TIR$condition=='pMM002.4 x CA1200 - 0.25% eth' | df_TIR$condition=='pMM002.4 x CA1200 - 250 uM auxin '),3:10], function(i) {wilcox.test(i ~ df_TIR$condition[which(df_TIR$condition=='pMM002.4 x CA1200 - 0.25% eth' | df_TIR$condition=='pMM002.4 x CA1200 - 250 uM auxin ')])$p.value})
pval_63uM <- sapply(df_TIR[which(df_TIR$condition=='pMM002.4 x CA1200 - 0.25% eth' | df_TIR$condition=='pMM002.4 x CA1200 - 62.5 uM auxin '),3:10], function(i) {wilcox.test(i ~ df_TIR$condition[which(df_TIR$condition=='pMM002.4 x CA1200 - 0.25% eth' | df_TIR$condition=='pMM002.4 x CA1200 - 62.5 uM auxin ')])$p.value})
pval_16uM <- sapply(df_TIR[which(df_TIR$condition=='pMM002.4 x CA1200 - 0.25% eth' | df_TIR$condition=='pMM002.4 x CA1200 - 15.6 uM auxin '),3:10], function(i) {wilcox.test(i ~ df_TIR$condition[which(df_TIR$condition=='pMM002.4 x CA1200 - 0.25% eth' | df_TIR$condition=='pMM002.4 x CA1200 - 15.6 uM auxin ')])$p.value})
pval_4uM <- sapply(df_TIR[which(df_TIR$condition=='pMM002.4 x CA1200 - 0.25% eth' | df_TIR$condition=='pMM002.4 x CA1200 - 3.9 uM auxin '),3:10], function(i) {wilcox.test(i ~ df_TIR$condition[which(df_TIR$condition=='pMM002.4 x CA1200 - 0.25% eth' | df_TIR$condition=='pMM002.4 x CA1200 - 3.9 uM auxin ')])$p.value})

pval_TIR <- data.frame (pval_4mM=pval_4mM, pval_1mM=pval_1mM, pval_250uM=pval_250uM, pval_63uM=pval_63uM,pval_16uM=pval_16uM,pval_4uM=pval_4uM)
write.csv(pval_TIR,'Plots/2020-11-10 TIR1only_pval.csv')

