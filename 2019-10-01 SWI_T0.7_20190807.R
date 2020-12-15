writeLines(capture.output(sessionInfo()), "2019- _sessionInfo.txt")

library(reshape2)
library(graphics)
library(scales)

setwd("~/Milou - Projects/GRH-1/SWI/2019-08-07")
Molt <- read.csv('2019-08-07 MoltAnnotation.csv')
GFP <- read.csv('Results T07/2019-08-07 GFP intensities T07.csv')
GFP_BGsub <- GFP[,-3] #take only background subtracted GFP intensities
tp <- 98
strain1 <- 'grh-1'
strain2 <- 'qua-1'

#data formatting
infoMolt <- Molt[Molt$valid==1,] #only valid positions
infoMolt$M3entry_M2exit <- infoMolt$M3entry - infoMolt$M2exit + 1 #correct time of M3 entry for M2 exit
infoMolt$M3exit_M2exit <- infoMolt$M3exit - infoMolt$M2exit + 1 #correct time of M3 exit for M2 exit
infoMolt$phenotype_M2exit <- infoMolt$phenotype - infoMolt$M2exit + 1 #correct time of phenotype for M2 exit

gfpInten <- dcast(GFP_BGsub,Position~Frame) #transform into matrix format
gfpInten <- gfpInten[gfpInten$Position %in% infoMolt$position,] #take only valid positions
rownames(gfpInten) <- gfpInten$Position
gfpInten <- gfpInten[,-1] #remove position 
colnames(gfpInten) <- (as.numeric(colnames(gfpInten))/6) - (1/6)

#function to plot GFP intensity relative to M2 exit (grey), mean GFP(black), SEM(red,dashed), M3(blue dashed)
plotGFPmolt <- function(gfp, molt,tp, name){
  #align GFP intensities by M2exit
  gfp_M2exit <- matrix(NA,nrow=nrow(gfp), ncol=tp)
  for (i in 1:nrow(gfp)){
    start <- molt$M2exit[i]
    end <- tp 
    gfp_M2exit[i,1:length(start:end)] <- as.matrix(gfp[i,start:end])
  }
  
  #replace data with NA if nr of worms smaller or equal to 5
  idx <- apply(gfp_M2exit,2,function(x){sum(!is.na(x))})
  gfp_M2exit[,which(idx<=5)] <-NA
  
  #replace data with NA for mean+SEM GFP
  ind <- which(gfp_M2exit < 50, arr.ind = T)
  gfp_M2exit_c <- gfp_M2exit
  gfp_M2exit_c[ind] <- NA
  
  #calculate mean+SEM of GFP intensity over all worms for each timepoint
  meanGFP <- c()
  sdGFP <- c()
  semGFP <- c()
  for (i in 1:ncol(gfp_M2exit_c)){
    meanGFP[i] <- mean(gfp_M2exit_c[,i], na.rm=TRUE)
    sdGFP[i] <- sd(gfp_M2exit_c[,i], na.rm=TRUE)
    semGFP[i] <- sd(gfp_M2exit_c[,i], na.rm=TRUE)/sqrt(sum(!is.na(gfp_M2exit_c[,i])))
  }
  
  upperGFP <- meanGFP+semGFP
  lowerGFP <- meanGFP-semGFP
  
  #correct M3 entry and M3 exit for M2 exit
  gfp_molt <- matrix(NA,nrow=nrow(molt), ncol=tp)
  for (i in 1:nrow(molt)){
    gfp_molt[i,molt$M3entry_M2exit[i]:molt$M3exit_M2exit[i]] <- gfp_M2exit[i,molt$M3entry_M2exit[i]:molt$M3exit_M2exit[i]]
  }
  
  return(list(meanGFP=meanGFP,upperGFP=upperGFP,lowerGFP=lowerGFP, gfp_M2exit=gfp_M2exit, gfp_molt=gfp_molt, name=name, gfp=gfp))
}

plotGFPmolt_aux <- function(gfp, molt,tp, name){
  #align GFP intensities by M2exit
  gfp_M2exit <- matrix(NA,nrow=nrow(gfp), ncol=tp)
  for (i in 1:nrow(gfp)){
    start <- molt$M2exit[i]
    end <- tp 
    gfp_M2exit[i,1:length(start:end)] <- as.matrix(gfp[i,start:end])
  }
  
  #replace data with NA after bursting
  for (i in 1:nrow(molt)){
    gfp_M2exit[i,molt$phenotype_M2exit[i]:tp] <- NA
  }
  
  #replace data with NA if nr of worms smaller or equal to 5
  idx <- apply(gfp_M2exit,2,function(x){sum(!is.na(x))})
  gfp_M2exit[,which(idx<=5)] <-NA
  
  #replace data with NA for mean+SEM GFP
  ind <- which(gfp_M2exit < 50, arr.ind = T)
  gfp_M2exit_c <- gfp_M2exit
  gfp_M2exit_c[ind] <- NA
  
  #calculate mean+SEM of GFP intensity over all worms for each timepoint
  meanGFP <- c()
  sdGFP <- c()
  semGFP <- c()
  for (i in 1:ncol(gfp_M2exit_c)){
    meanGFP[i] <- mean(gfp_M2exit_c[,i], na.rm=TRUE)
    sdGFP[i] <- sd(gfp_M2exit_c[,i], na.rm=TRUE)
    semGFP[i] <- sd(gfp_M2exit_c[,i], na.rm=TRUE)/sqrt(sum(!is.na(gfp_M2exit_c[,i])))
  }
  
  upperGFP <- meanGFP+semGFP
  lowerGFP <- meanGFP-semGFP
  
  #correct M3 entry and phenotype for M2 exit
  gfp_molt_phenotype <- matrix(NA,nrow=nrow(molt), ncol=tp)
  for (i in 1:nrow(molt)){
    gfp_molt_phenotype[i,molt$M3entry_M2exit[i]:molt$phenotype_M2exit[i]] <- gfp_M2exit[i,molt$M3entry_M2exit[i]:molt$phenotype_M2exit[i]]
  }
  
  return(list(meanGFP=meanGFP,upperGFP=upperGFP,lowerGFP=lowerGFP, gfp_M2exit=gfp_M2exit, gfp_molt_phenotype=gfp_molt_phenotype, name=name, gfp=gfp))
}


GFP_strain1_eth <- plotGFPmolt(gfpInten[infoMolt$condition=='ethanol' & infoMolt$strain==strain1,],infoMolt[infoMolt$condition=='ethanol' & infoMolt$strain==strain1,],tp, paste(strain1, 'ethanol'))
GFP_strain2_eth <- plotGFPmolt(gfpInten[infoMolt$condition=='ethanol' & infoMolt$strain==strain2,],infoMolt[infoMolt$condition=='ethanol' & infoMolt$strain==strain2,],tp, paste(strain2, 'ethanol'))

GFP_strain1_aux <- plotGFPmolt_aux(gfpInten[infoMolt$condition=='auxin' & infoMolt$strain==strain1 & infoMolt$phen.descrip=='L3 phenotype',],infoMolt[infoMolt$condition=='auxin' & infoMolt$strain==strain1 & infoMolt$phen.descrip=='L3 phenotype',],tp, paste(strain1, 'auxin'))
GFP_strain2_aux <- plotGFPmolt_aux(gfpInten[infoMolt$condition=='auxin' & infoMolt$strain==strain2 & infoMolt$phen.descrip=='L3 phenotype',],infoMolt[infoMolt$condition=='auxin' & infoMolt$strain==strain2 & infoMolt$phen.descrip=='L3 phenotype',],tp, paste(strain2, 'auxin'))

plot(GFP_strain1_eth$meanGFP, type='l', main=strain1, ylim=c(0,1100), xaxt='n'); lines(GFP_strain1_aux$meanGFP, col='red'); lines(GFP_strain1_eth$upperGFP, col='grey'); lines(GFP_strain1_eth$lowerGFP, col='grey'); lines(GFP_strain1_aux$upperGFP, col='pink'); lines(GFP_strain1_aux$lowerGFP, col='pink'); axis(1, at = 1:tp, labels = colnames(GFP_strain1_eth$gfp), cex.axis = 0.7); abline(v=mean(infoMolt$M3entry_M2exit[which(infoMolt$condition=='ethanol' & infoMolt$strain=='grh-1')])); abline(v=mean(infoMolt$M3exit_M2exit[which(infoMolt$condition=='ethanol' & infoMolt$strain=='grh-1')])); abline(v=mean(infoMolt$M3entry_M2exit[which(infoMolt$condition=='auxin' & infoMolt$strain=='grh-1' & infoMolt$phen.descrip=='L3 phenotype')]), col='red')
plot(GFP_strain2_eth$meanGFP, type='l', main=strain2, ylim=c(0,500), xaxt='n'); lines(GFP_strain2_aux$meanGFP, col='red'); lines(GFP_strain2_eth$upperGFP, col='grey'); lines(GFP_strain2_eth$lowerGFP, col='grey'); lines(GFP_strain2_aux$upperGFP, col='pink'); lines(GFP_strain2_aux$lowerGFP, col='pink'); axis(1, at = 1:tp, labels = colnames(GFP_strain1_eth$gfp), cex.axis = 0.7); abline(v=mean(infoMolt$M3entry_M2exit[which(infoMolt$condition=='ethanol' & infoMolt$strain=='qua-1')])); abline(v=mean(infoMolt$M3exit_M2exit[which(infoMolt$condition=='ethanol' & infoMolt$strain=='qua-1')])); abline(v=mean(infoMolt$M3entry_M2exit[which(infoMolt$condition=='auxin' & infoMolt$strain=='qua-1' & infoMolt$phen.descrip=='L3 phenotype')]), col='red')

matplot(t(GFP_strain1_eth$gfp_M2exit), ylim=c(0,1100), type='l', lty=1, col='grey', main=GFP_strain1_eth$name,xaxt = 'n', ylab='background subtracted GFP intensity', xlab='time after M2 exit (h)'); matplot(t(GFP_strain1_eth$gfp_molt), type='l', lty=2, col=alpha('blue',0.5), add=TRUE); lines(GFP_strain1_eth$meanGFP);  lines(GFP_strain1_eth$upperGFP, col='red', lty=2); lines(GFP_strain1_eth$lowerGFP, col='red', lty=2); axis(1, at = 1:tp, labels = colnames(GFP_strain1_eth$gfp), cex.axis = 0.7) 
matplot(t(GFP_strain1_aux$gfp_M2exit), ylim=c(0,1100), type='l', lty=1, col='grey', main=GFP_strain1_aux$name,xaxt = 'n', ylab='background subtracted GFP intensity', xlab='time after M2 exit (h)'); matplot(t(GFP_strain1_aux$gfp_molt_phenotype), type='l', lty=2, col=alpha('blue',0.5), add=TRUE); lines(GFP_strain1_aux$meanGFP);  lines(GFP_strain1_aux$upperGFP, col='red', lty=2); lines(GFP_strain1_aux$lowerGFP, col='red', lty=2); axis(1, at = 1:tp, labels = colnames(GFP_strain1_aux$gfp), cex.axis = 0.7) 
matplot(t(GFP_strain2_eth$gfp_M2exit), ylim=c(0,500), type='l', lty=1, col='grey', main=GFP_strain2_eth$name,xaxt = 'n', ylab='background subtracted GFP intensity', xlab='time after M2 exit (h)'); matplot(t(GFP_strain2_eth$gfp_molt), type='l', lty=2, col=alpha('blue',0.5), add=TRUE); lines(GFP_strain2_eth$meanGFP);  lines(GFP_strain2_eth$upperGFP, col='red', lty=2); lines(GFP_strain2_eth$lowerGFP, col='red', lty=2); axis(1, at = 1:tp, labels = colnames(GFP_strain2_eth$gfp), cex.axis = 0.7) 
matplot(t(GFP_strain2_aux$gfp_M2exit), ylim=c(0,500), type='l', lty=1, col='grey', main=GFP_strain2_aux$name,xaxt = 'n', ylab='background subtracted GFP intensity', xlab='time after M2 exit (h)'); matplot(t(GFP_strain2_aux$gfp_molt_phenotype), type='l', lty=2, col=alpha('blue',0.5), add=TRUE); lines(GFP_strain2_aux$meanGFP);  lines(GFP_strain2_aux$upperGFP, col='red', lty=2); lines(GFP_strain2_aux$lowerGFP, col='red', lty=2); axis(1, at = 1:tp, labels = colnames(GFP_strain2_aux$gfp), cex.axis = 0.7) 

#mean and SEM do not include measurements for which background subtracted GFP intensity < 50
#plot is cut when nr of worms is smaller or equal to 5 
mean(infoMolt$M3entry_M2exit[which(infoMolt$condition=='ethanol' & infoMolt$strain=='grh-1')])/6
#8.231884
mean(infoMolt$M3exit_M2exit[which(infoMolt$condition=='ethanol' & infoMolt$strain=='grh-1')])/6
#10.65942
mean(infoMolt$M3entry_M2exit[which(infoMolt$condition=='ethanol' & infoMolt$strain=='qua-1')])/6
#7.653846
mean(infoMolt$M3exit_M2exit[which(infoMolt$condition=='ethanol' & infoMolt$strain=='qua-1')])/6
#10.0641

