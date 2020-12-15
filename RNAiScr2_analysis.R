writeLines(capture.output(sessionInfo()), "RNAiScr2_analysis_sessionInfo.txt")

library('R.matlab')
library('reshape2')
library('ggplot2')
library('survminer')
library('stats')
library("RColorBrewer")
setwd("/tungstenfs/scratch/ggrossha/meeumilo/scripts/ggrosshans_grh-1_paper")
set.seed(123)

#make barplot of completion of molts
moltCompl <- read.table('Data/2019-03-25 quantification_moltcompletion_hits.csv', sep=',',skip=2, header=TRUE, row.names=1)
colnames(moltCompl) <- c('M1', 'M1-M2', 'M1-M3', 'M1-M4', 'total')
moltComplPer <- moltCompl[,1:4]/moltCompl$total
barplot(as.matrix(t(moltComplPer[,2:4])), legend = colnames(moltComplPer[,2:4]), main='Completion of Molts', ylab='% of animals', xlab='RNAi')

#quantify and plot duration of stages
file_Scr01 <- readMat('Data/2017-06-12 RNAiScr01.mat')
file_Scr02 <- readMat('Data/2017-06-16 RNAiScr02.mat')
file_Scr08 <- readMat('Data/2017-06-26 RNAiScr08.mat')

myfunction <- function(file){
  
  valid <- as.factor(file$valid)
  annotation <- unlist(file$samples)
  
  M1 <- c()
  M2 <- c()
  M3 <- c()
  M4 <- c()
  I1 <- c()
  I2 <- c()
  I3 <- c()
  I4 <- c()
  
  for (i in 1:ncol(file$Hatch)){
    M1[i] <- file$Molt[1,2,i] - file$Molt[1,1,i]
    M2[i] <- file$Molt[2,2,i] - file$Molt[2,1,i]
    M3[i] <- file$Molt[3,2,i] - file$Molt[3,1,i]
    M4[i] <- file$Molt[4,2,i] - file$Molt[4,1,i]
    I1[i] <- file$Molt[1,1,i] - file$Hatch[1,i]
    I2[i] <- file$Molt[2,1,i] - file$Molt[1,2,i]
    I3[i] <- file$Molt[3,1,i] - file$Molt[2,2,i]
    I4[i] <- file$Molt[4,1,i] - file$Molt[3,2,i]
  }
  
  df <- data.frame(condition = annotation, valid = valid, I1=I1,M1=M1,I2=I2,M2=M2,I3=I3,M3=M3,I4=I4, M4=M4)
  df <- df[which(df$valid==1),]
  return(df)
  }

df_Scr01 <- myfunction(file_Scr01)
df_Scr02 <- myfunction(file_Scr02)
df_Scr08 <- myfunction(file_Scr08)


df_all <- rbind(df_Scr01,df_Scr02,df_Scr08)
df_all$condition <- ordered(df_all$condition, levels=c('EV - rde', 'EV', 'bed3 - rde', 'bed3', 'blmp - rde', 'blmp', 'nhr25 - rde', 'nhr25')) #order the annotation to get right order in boxplot


#plots
par(mfrow=c(2,2))
boxplot(df_all$M1/6 ~ df_all$condition, data=df_all, ylim=c(0,7.5), main='Molt 1', las=2, col=rep(c('azure2','azure4'),4))
boxplot(df_all$M2/6 ~ df_all$condition, data=df_all, ylim=c(0,7.5), main='Molt 2', las=2, col=rep(c('azure2','azure4'),4))
boxplot(df_all$M3/6 ~ df_all$condition, data=df_all, ylim=c(0,7.5), main='Molt 3', las=2, col=rep(c('azure2','azure4'),4))
boxplot(df_all$M4/6 ~ df_all$condition, data=df_all, ylim=c(0,7.5), main='Molt 4', las=2, col=rep(c('azure2','azure4'),4))
boxplot(df_all$I1/6 ~ df_all$condition, data=df_all, ylim=c(4,16), main='Intermolt 1', las=2, col=rep(c('azure2','azure4'),4))
boxplot(df_all$I2/6 ~ df_all$condition, data=df_all, ylim=c(4,16), main='Intermolt 2', las=2, col=rep(c('azure2','azure4'),4))
boxplot(df_all$I3/6 ~ df_all$condition, data=df_all, ylim=c(4,16), main='Intermolt 3', las=2, col=rep(c('azure2','azure4'),4))
boxplot(df_all$I4/6 ~ df_all$condition, data=df_all, ylim=c(4,16), main='Intermolt 4', las=2, col=rep(c('azure2','azure4'),4))

#pval
df_Scr02_EV <- df_Scr02[grepl('EV', df_Scr02$condition),]
df_Scr02_blmp <- df_Scr02[grepl('blmp', df_Scr02$condition),]
wilcox.test(df_all$I4[which(df_all$condition=='bed3')], df_all$I4[which(df_all$condition=='bed3 - rde')])$p.value
pval_bed3 <- sapply(df_Scr01[,3:10], function(i) {wilcox.test(i ~ df_Scr01$condition)$p.value})
pval_EV <- sapply(df_Scr02_EV[,3:10], function(i) {wilcox.test(i ~ df_Scr02_EV$condition)$p.value})
pval_blmp <- sapply(df_Scr02_blmp[,3:10], function(i) {wilcox.test(i ~ df_Scr02_blmp$condition)$p.value})
pval_nhr25 <- sapply(df_Scr08[,3:10], function(i) {wilcox.test(i ~ df_Scr08$condition)$p.value})

pval <- data.frame(EV = pval_EV, bed3 = pval_bed3, blmp = pval_blmp, nhr25 = pval_nhr25) 
write.csv(pval,'Plots/LucScreen/2020-11-10 pval_boxplots_molt_intermolt.csv')

#fold change
FC_bed3 <- (colMeans(df_Scr01[which(df_Scr01$condition=='bed3'),3:10])/6) / (colMeans(df_Scr01[which(df_Scr01$condition=='bed3 - rde'),3:10])/6)
FC_nhr25 <- (colMeans(df_Scr08[which(df_Scr08$condition=='nhr25'),3:10])/6) / (colMeans(df_Scr08[which(df_Scr08$condition=='nhr25 - rde'),3:10])/6)
FC_EV <- (colMeans(df_Scr02_EV[which(df_Scr02_EV$condition=='EV'),3:10])/6) / (colMeans(df_Scr02_EV[which(df_Scr02_EV$condition=='EV - rde'),3:10])/6)
FC_blmp <- (colMeans(df_Scr02_blmp[which(df_Scr02_blmp$condition=='blmp'),3:10])/6) / (colMeans(df_Scr02_blmp[which(df_Scr02_blmp$condition=='blmp - rde'),3:10])/6)

FC <- data.frame(EV=FC_EV, bed3=FC_bed3,blmp=FC_blmp,nhr25=FC_nhr25)







