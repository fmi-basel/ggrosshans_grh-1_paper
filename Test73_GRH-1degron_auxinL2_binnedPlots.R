writeLines(capture.output(sessionInfo()), "Test73_binnedPlots_sessionInfo.txt")

'
Analyze Luciferase Assay data Test73

DESCRIPTION
-----------
This code allows the user to analyze the luciferase assay data of Test73 and generate figures of the GRH-1 manuscript

This script consists of the following steps: 
  (1) load data
  (2) calculate stage durations
  (3) group according to annotated phenotypes
  (4) make plots: scatter, box, waterfall, heatmap


REFERENCES
----------
LuciferaseAssayAnalyzer (Meeuse, 2020)


USAGE
-----
Data format: 
file: .mat file from LuciferaseAssayAnalyzer containing molt and hatch annotations, luminescence values, and sample names
Phenot_annot: .csv file containing mannual annotation of phenotype, i.e. whether in this case scored 0 when M2 was last observed molt and scored 1 when M3 was last observed molt


INPUT
-----
2018-12-03 auxin in L2.mat
2019-08-27 Test73_QuantificationM2execution.csv

OUTPUT
------
plots: scatter, box, waterfall, heatmap


REQUIREMENTS
------------
This script requires the following packages to be installed:
R version 4.0.0 (2020-04-24)
cowplot_1.0.0  ggplot2_3.3.0  R.matlab_3.6.2


"'

library('R.matlab')
library(ggplot2)
library(cowplot)
library(NMF)
library('zoo')

setwd("/tungstenfs/scratch/ggrossha/meeumilo/scripts/ggrosshans_grh-1_paper")

#---Load data---
file <- readMat('Data/2018-12-03 auxin in L2.mat')
valid <- as.factor(file$valid)
annotation <- unlist(file$samples)

#---calculate stage durations---
M1 <- c()
M2 <- c()
I1 <- c()
I2 <- c()
L1 <- c()
L2 <- c()
M1_in <- c()
M2_in <- c()
M3_in <- c()
M4_in <- c()
M3 <- c()
M4 <- c()
I3 <- c()
I4 <- c()
L3 <- c()
L4 <- c()
hatch <- c()

for (i in 1:length(valid)){
  M1[i] <- (file$Molt[1,2,i] - file$Molt[1,1,i])/6
  M2[i] <- (file$Molt[2,2,i] - file$Molt[2,1,i])/6
  I1[i] <- (file$Molt[1,1,i] - file$Hatch[1,i])/6
  I2[i] <- (file$Molt[2,1,i] - file$Molt[1,2,i])/6
  L1[i] <- (file$Molt[1,2,i] - file$Hatch[1,i])/6
  L2[i] <- (file$Molt[2,2,i] - file$Molt[1,2,i])/6
  M1_in[i] <- file$Molt[1,1,i]/6 #time of M1 entry
  M2_in[i] <- file$Molt[2,1,i]/6
  M3_in[i] <- file$Molt[3,1,i]/6
  M4_in[i] <- file$Molt[4,1,i]/6
  M3[i] <- (file$Molt[3,2,i] - file$Molt[3,1,i])/6
  M4[i] <- (file$Molt[4,2,i] - file$Molt[4,1,i])/6
  I3[i] <- (file$Molt[3,1,i] - file$Molt[2,2,i])/6
  I4[i] <- (file$Molt[4,1,i] - file$Molt[3,2,i])/6
  L3[i] <- (file$Molt[3,2,i] - file$Molt[2,2,i])/6
  L4[i] <- (file$Molt[4,2,i] - file$Molt[3,2,i])/6
  hatch[i] <- file$Hatch[1,i]/6
}


df <- data.frame(condition = annotation, valid = valid, hatch=hatch, I1=I1,M1=M1, L1=L1, M1_in=M1_in, I2=I2,M2=M2, L2=L2, M2_in=M2_in,I3=I3,M3=M3, L3=L3, M3_in=M3_in, I4=I4,M4=M4, L4=L4, M4_in=M4_in)
df <- df[which(valid==1),] #take only samples annotated as valid
df_aux <- subset(df, condition=='auxin')


#---group according to phenotype in M2 ('0') vs M3 ('1')---
Phenot_annot <- read.csv('Data/2019-08-27 Test73_QuantificationM2execution.csv', row.names = 1)
Phenot_annot <- Phenot_annot[which(valid==1),]
Phenot_annot <- as.factor(Phenot_annot)
df$phenAnnot <- Phenot_annot #add phenotype annotation to df

df_M2phenot <- df[which(Phenot_annot==0),]
df_M3phenot <- df[which(Phenot_annot==1),] #controls are also annotated with 1, subset further to extract aux
df_M2phenot_aux <- subset(df_M2phenot, condition=='auxin')
df_M3phenot_aux <- subset(df_M3phenot, condition=='auxin')
df_ctrl <- subset(df, condition=='vehicle')
df_ctrl_sel <- df_ctrl[which(df_ctrl$M2_in<34),]

#---scatterplot of time of M1 entry vs stage durations for M2 phenotype---
plotM2_phenot_I1 <- ggplot(df_M2phenot_aux) + geom_point(aes(x=M1_in, y=I1), alpha=0.25) + labs(x='Time of M1 entry (h)', y='Duration of I1 (h)', title='phenotype in M2')+ylim(min(df$I1)*0.95,max(df$I1)*1.025) +xlim(min(df$M1_in)*0.95,max(df$M1_in)*1.025) +geom_vline(xintercept=29) + theme_classic()
plotM2_phenot_M1 <- ggplot(df_M2phenot_aux) + geom_point(aes(x=M1_in, y=M1), alpha=0.25) + labs(x='Time of M1 entry (h)', y='Duration of M1 (h)', title='phenotype in M2')+ylim(min(df$M1)*0.95,max(df$M1)*1.025) +xlim(min(df$M1_in)*0.95,max(df$M1_in)*1.025) +geom_vline(xintercept=29) + theme_classic()
plotM2_phenot_I2 <- ggplot(df_M2phenot_aux) + geom_point(aes(x=M1_in, y=I2), alpha=0.25) + labs(x='Time of M1 entry (h)', y='Duration of I2 (h)', title='phenotype in M2')+ylim(min(df$I2)*0.95,max(df$I2)*1.025) +xlim(min(df$M1_in)*0.95,max(df$M1_in)*1.025) +geom_vline(xintercept=29) + theme_classic()
plotM3_phenot_ctrl_I1 <- ggplot(df_ctrl) + geom_point(aes(x=M1_in, y=I1), alpha=0.25) + labs(x='Time of M1 entry (h)', y='Duration of I1 (h)', title='control')+ylim(min(df$I1)*0.95,max(df$I1)*1.025) +xlim(min(df$M1_in)*0.95,max(df$M1_in)*1.025) +geom_vline(xintercept=29) + theme_classic()
plotM3_phenot_ctrl_M1 <- ggplot(df_ctrl) + geom_point(aes(x=M1_in, y=M1), alpha=0.25) + labs(x='Time of M1 entry (h)', y='Duration of M1 (h)', title='control')+ylim(min(df$M1)*0.95,max(df$M1)*1.025) +xlim(min(df$M1_in)*0.95,max(df$M1_in)*1.025) +geom_vline(xintercept=29) + theme_classic()
plotM3_phenot_ctrl_I2 <- ggplot(df_ctrl) + geom_point(aes(x=M1_in, y=I2), alpha=0.25) + labs(x='Time of M1 entry (h)', y='Duration of I2 (h)', title='control')+ylim(min(df$I2)*0.95,max(df$I2)*1.025) +xlim(min(df$M1_in)*0.95,max(df$M1_in)*1.025) +geom_vline(xintercept=29) + theme_classic()
plot_grid(plotM2_phenot_I1, plotM3_phenot_ctrl_I1,plotM2_phenot_M1, plotM3_phenot_ctrl_M1,plotM2_phenot_I2, plotM3_phenot_ctrl_I2, nrow=3, ncol=2, labels = 'AUTO')

#---scatterplot of time of M2 entry vs stage durations for M3 phenotype---
plotM3_phenot_M2 <- ggplot(df_M3phenot_aux) + geom_point(aes(x=M2_in, y=M2), alpha=0.25) + labs(x='Time of M2 entry (h)', y='Duration of M2 (h)', title='phenotype in M3')+ylim(min(c(df_ctrl$M2, df_M3phenot_aux$M2))*0.95,max(c(df_ctrl$M2, df_M3phenot_aux$M2))*1.025) +xlim(min(df$M2_in)*0.95,max(df$M2_in)*1.025) +geom_vline(xintercept=29) + theme_classic()
plotM3_phenot_I2 <- ggplot(df_M3phenot_aux) + geom_point(aes(x=M2_in, y=I2), alpha=0.25) + labs(x='Time of M2 entry (h)', y='Duration of I2 (h)', title='phenotype in M3')+ylim(min(c(df_ctrl$I2, df_M3phenot_aux$I2))*0.95,max(c(df_ctrl$I2, df_M3phenot_aux$I2))*1.025) +xlim(min(df$M2_in)*0.95,max(df$M2_in)*1.025) +geom_vline(xintercept=29) + theme_classic()
plotM3_phenot_I3 <- ggplot(df_M3phenot_aux) + geom_point(aes(x=M2_in, y=I3), alpha=0.25) + labs(x='Time of M2 entry (h)', y='Duration of I3 (h)', title='phenotype in M3')+ylim(min(c(df_ctrl$I3, df_M3phenot_aux$I3))*0.95,max(c(df_ctrl$I3, df_M3phenot_aux$I3))*1.025) +xlim(min(df$M2_in)*0.95,max(df$M2_in)*1.025) +geom_vline(xintercept=29) + theme_classic()
plotM3_phenot_ctrl_M2 <- ggplot(df_ctrl) + geom_point(aes(x=M2_in, y=M2), alpha=0.25) + labs(x='Time of M2 entry (h)', y='Duration of M2 (h)', title='control')+ylim(min(c(df_ctrl$M2, df_M3phenot_aux$M2))*0.95,max(c(df_ctrl$M2, df_M3phenot_aux$M2))*1.025) +xlim(min(df$M2_in)*0.95,max(df$M2_in)*1.025) +geom_vline(xintercept=29) + theme_classic()
plotM3_phenot_ctrl_I2 <- ggplot(df_ctrl) + geom_point(aes(x=M2_in, y=I2), alpha=0.25) + labs(x='Time of M2 entry (h)', y='Duration of I2 (h)', title='control')+ylim(min(c(df_ctrl$I2, df_M3phenot_aux$I2))*0.95,max(c(df_ctrl$I2, df_M3phenot_aux$I2))*1.025) +xlim(min(df$M2_in)*0.95,max(df$M2_in)*1.025) +geom_vline(xintercept=29) + theme_classic()
plotM3_phenot_ctrl_I3 <- ggplot(df_ctrl) + geom_point(aes(x=M2_in, y=I3), alpha=0.25) + labs(x='Time of M2 entry (h)', y='Duration of I3 (h)', title='control')+ylim(min(c(df_ctrl$I3, df_M3phenot_aux$I3))*0.95,max(c(df_ctrl$I3, df_M3phenot_aux$I3))*1.025) +xlim(min(df$M2_in)*0.95,max(df$M2_in)*1.025) +geom_vline(xintercept=29) + theme_classic()
plot_grid(plotM3_phenot_I2, plotM3_phenot_ctrl_I2,plotM3_phenot_M2, plotM3_phenot_ctrl_M2, plotM3_phenot_I3, plotM3_phenot_ctrl_I3,nrow=3, ncol=2, labels = 'AUTO')


#---scatterplot with binned mean and stdev---
functionDurBin <- function(df.toMean, df.toBin, df.toBreaks){
  binBreaks <- seq(floor(min(df.toBreaks)),ceiling(max(df.toBreaks)),by=1)
  binned <-.bincode(df.toBin, breaks=binBreaks, right=TRUE, include.lowest=TRUE) #right: intervals should be closed on the right and open on the left 
  meanDur.bin <- data.frame(bin=binBreaks) 
  for (i in 1:length(binBreaks)){
    meanDur.bin$stdev[i] <- sd(df.toMean[binned==i])
    meanDur.bin$mean[i] <- mean(df.toMean[binned==i])
    meanDur.bin$mean[is.nan(meanDur.bin$mean)] <- NA #set nan to NA
    meanDur.bin$mean[is.na(meanDur.bin$stdev)] <- NA #exclude mean of bins with only 1 datapoint -> the mean and stdev for this bin will not be plotted, but the single datapoint will be
  }
  meanDur.bin <- meanDur.bin[!is.na(meanDur.bin$mean),] #remove bins with NA
  return(meanDur.bin)
  }

M3phen_M2binned <- functionDurBin(df_M3phenot_aux$M2, df_M3phenot_aux$M2_in, c(df_M3phenot_aux$M2_in,df_ctrl_sel$M2_in))
M3phen_I2binned <- functionDurBin(df_M3phenot_aux$I2, df_M3phenot_aux$M2_in, c(df_M3phenot_aux$M2_in,df_ctrl_sel$M2_in))
M3phen_I3binned <- functionDurBin(df_M3phenot_aux$I3, df_M3phenot_aux$M2_in, c(df_M3phenot_aux$M2_in,df_ctrl_sel$M2_in))
M3ctl_M2binned <- functionDurBin(df_ctrl_sel$M2, df_ctrl_sel$M2_in, c(df_M3phenot_aux$M2_in,df_ctrl_sel$M2_in))
M3ctl_I2binned <- functionDurBin(df_ctrl_sel$I2, df_ctrl_sel$M2_in, c(df_M3phenot_aux$M2_in,df_ctrl_sel$M2_in))
M3ctl_I3binned <- functionDurBin(df_ctrl_sel$I3, df_ctrl_sel$M2_in, c(df_M3phenot_aux$M2_in,df_ctrl_sel$M2_in))

#plots: mean and stdev of bins with 1 datapoint are not shown
par(mfrow=c(2,3))

xlim_plot <- c(0.95*min(c(df_M3phenot_aux$M2_in,df_ctrl_sel$M2_in)),1.05*max(c(df_M3phenot_aux$M2_in,df_ctrl_sel$M2_in)))

plot(df_M3phenot_aux$M2_in, df_M3phenot_aux$I2, col=alpha('turquoise',0.2), pch=20, xlim=xlim_plot); 
lines(M3phen_I2binned$bin+0.5, M3phen_I2binned$mean, col='turquoise', lwd=2);
points(df_ctrl_sel$M2_in, df_ctrl_sel$I2, col=alpha('black',0.2), pch=20); 
lines(M3ctl_I2binned$bin+0.5, M3ctl_I2binned$mean, col='black', lwd=2)
polygon(c(M3phen_I2binned$bin+0.5,rev(M3phen_I2binned$bin+0.5)), c(M3phen_I2binned$mean+M3phen_I2binned$stdev,rev(M3phen_I2binned$mean-M3phen_I2binned$stdev)), col=alpha('turquoise',0.05), border=NA)
polygon(c(M3ctl_I2binned$bin+0.5,rev(M3ctl_I2binned$bin+0.5)), c(M3ctl_I2binned$mean+M3ctl_I2binned$stdev,rev(M3ctl_I2binned$mean-M3ctl_I2binned$stdev)), col=alpha('black',0.05), border=NA)
#lines(M3ctl_I2binned$bin+0.5,M3ctl_I2binned$mean-M3ctl_I2binned$stdev, lty=2);
#lines(M3ctl_I2binned$bin+0.5,M3ctl_I2binned$mean+M3ctl_I2binned$stdev, lty=2);
#lines(M3phen_I2binned$bin+0.5,M3phen_I2binned$mean-M3phen_I2binned$stdev, lty=2, col='turquoise');
#lines(M3phen_I2binned$bin+0.5,M3phen_I2binned$mean+M3phen_I2binned$stdev, lty=2, col='turquoise');

plot(df_M3phenot_aux$M2_in, df_M3phenot_aux$M2, col=alpha('turquoise',0.2), pch=20, xlim=xlim_plot); 
lines(M3phen_M2binned$bin+0.5, M3phen_M2binned$mean, col='turquoise', lwd=2);
points(df_ctrl_sel$M2_in, df_ctrl_sel$M2, col=alpha('black',0.2), pch=20); 
lines(M3ctl_M2binned$bin+0.5, M3ctl_M2binned$mean, col='black', lwd=2);
polygon(c(M3phen_M2binned$bin+0.5,rev(M3phen_M2binned$bin+0.5)), c(M3phen_M2binned$mean+M3phen_M2binned$stdev,rev(M3phen_M2binned$mean-M3phen_M2binned$stdev)), col=alpha('turquoise',0.05), border=NA)
polygon(c(M3ctl_M2binned$bin+0.5,rev(M3ctl_M2binned$bin+0.5)), c(M3ctl_M2binned$mean+M3ctl_M2binned$stdev,rev(M3ctl_M2binned$mean-M3ctl_M2binned$stdev)), col=alpha('black',0.05), border=NA)
#lines(M3ctl_M2binned$bin+0.5,M3ctl_M2binned$mean-M3ctl_M2binned$stdev, lty=2);
#lines(M3ctl_M2binned$bin+0.5,M3ctl_M2binned$mean+M3ctl_M2binned$stdev, lty=2);
#lines(M3phen_M2binned$bin+0.5,M3phen_M2binned$mean-M3phen_M2binned$stdev, lty=2, col='turquoise');
#lines(M3phen_M2binned$bin+0.5,M3phen_M2binned$mean+M3phen_M2binned$stdev, lty=2, col='turquoise');

plot(df_M3phenot_aux$M2_in, df_M3phenot_aux$I3, col=alpha('turquoise',0.2), pch=20, xlim=xlim_plot); 
lines(M3phen_I3binned$bin+0.5, M3phen_I3binned$mean, col='turquoise', lwd=2);
points(df_ctrl_sel$M2_in, df_ctrl_sel$I3, col=alpha('black',0.2), pch=20); 
lines(M3ctl_I3binned$bin+0.5, M3ctl_I3binned$mean, col='black', lwd=2);
polygon(c(M3phen_I3binned$bin+0.5,rev(M3phen_I3binned$bin+0.5)), c(M3phen_I3binned$mean+M3phen_I3binned$stdev,rev(M3phen_I3binned$mean-M3phen_I3binned$stdev)), col=alpha('turquoise',0.05), border=NA)
polygon(c(M3ctl_I3binned$bin+0.5,rev(M3ctl_I3binned$bin+0.5)), c(M3ctl_I3binned$mean+M3ctl_I3binned$stdev,rev(M3ctl_I3binned$mean-M3ctl_I3binned$stdev)), col=alpha('black',0.05), border=NA)
#lines(M3ctl_I3binned$bin+0.5,M3ctl_I3binned$mean-M3ctl_I3binned$stdev, lty=2);
#lines(M3ctl_I3binned$bin+0.5,M3ctl_I3binned$mean+M3ctl_I3binned$stdev, lty=2);
#lines(M3phen_I3binned$bin+0.5,M3phen_I3binned$mean-M3phen_I3binned$stdev, lty=2, col='turquoise');
#lines(M3phen_I3binned$bin+0.5,M3phen_I3binned$mean+M3phen_I3binned$stdev, lty=2, col='turquoise');




#---waterfall plot of M2 entry colored by phenotype---
df_M2entry <- data.frame(M2entry=c(df_M2phenot_aux$M2_in,df_M3phenot_aux$M2_in), phenotype=c(rep('M2', time=length(df_M2phenot_aux$M2_in)),rep('M3', time=length(df_M3phenot_aux$M2_in))))
df_M2entry_ordered <- df_M2entry[order(df_M2entry$M2entry),-3]
df_M2entry_ordered$index <- rev(c(1:nrow(df_M2entry_ordered)))
ggplot(df_M2entry_ordered, aes(x=M2entry, y=index, col=phenotype)) + geom_point(alpha=0.5) + scale_color_manual(values=c("green", "blue")) + labs(x='Time of M2 entry') +geom_vline(xintercept=29)
plot_grid(ggplot(df_M2entry_ordered, aes(x=M2entry, y=index, col=phenotype)) + geom_point(alpha=0.5) + scale_color_manual(values=c("green", "blue")) + labs(x='Time of M2 entry') +geom_vline(xintercept=29),nrow=1, ncol=4, labels = 'AUTO')

#---make heatmaps of luminescence for L2 phenotype for roughly 20 animals---
L2phen_df_aux <- df_aux[which(df_aux$M2_in>33 & df_aux$M2_in<34),] #select M2 phenotype by manual inspection of time points in lumi heatmap
L2phen_df_aux <- L2phen_df_aux[order(L2phen_df_aux$M2_in),]
L2phen <- as.numeric(rownames(L2phen_df_aux))

myCol <- colorRampPalette(c('black','white'))(100)
aheatmap(t(file$Xc[,L2phen]), Rowv = NA, Colv = NA, breaks=seq(-1,4,length.out = 51), color = myCol, main='L2 phenotype')

