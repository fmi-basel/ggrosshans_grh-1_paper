writeLines(capture.output(sessionInfo()), "2019- _sessionInfo.txt")

library(reshape2)
library(graphics)
library(scales)

setwd("~/Milou - Projects/GRH-1/SWI/2019-08-23")
Molt <- read.csv('2019-08-23 moltAnnotation.csv')
GFP33 <- read.csv('2019-08-23 Intensity_worm33.csv', row.names = 1)
GFP34 <- read.csv('2019-08-23 Intensity_worm34.csv', row.names = 1)
GFP60 <- read.csv('2019-08-23 Intensity_worm60.csv', row.names = 1)
tp <- 361
strain1 <- 'HW2603'


#data formatting
infoMolt <- Molt[Molt$valid==1 & Molt$phenotype==strain1,] #only valid positions
gfpInten33 <- c(GFP33$max, rep(NA, 361-nrow(GFP33)))
gfpInten34 <- c(GFP34$max, rep(NA, 361-nrow(GFP34)))
gfpInten60 <- c(GFP60$max, rep(NA, 361-nrow(GFP60)))
names(gfpInten33) <- (as.numeric(rownames(GFP34))/6) - (1/6)
names(gfpInten34) <- (as.numeric(rownames(GFP34))/6) - (1/6)
names(gfpInten60) <- (as.numeric(rownames(GFP34))/6) - (1/6)

#plot seperately
par(mfrow=c(3,1))
plot(seq(0,60,(1/6)),gfpInten33); abline(v=infoMolt$hatch[1]/6); abline(v=infoMolt$M1.entry[1]/6); abline(v=infoMolt$M1.exit[1]/6); abline(v=infoMolt$M2.entry[1]/6); abline(v=infoMolt$M2.exit[1]/6); abline(v=infoMolt$M3.entry[1]/6); abline(v=infoMolt$M3.exit[1]/6); abline(v=infoMolt$M4.entry[1]/6); abline(v=infoMolt$M4.exit[1]/6)
plot(seq(0,60,(1/6)),gfpInten34); abline(v=infoMolt$hatch[2]/6); abline(v=infoMolt$M1.entry[2]/6); abline(v=infoMolt$M1.exit[2]/6); abline(v=infoMolt$M2.entry[2]/6); abline(v=infoMolt$M2.exit[2]/6); abline(v=infoMolt$M3.entry[2]/6); abline(v=infoMolt$M3.exit[2]/6); abline(v=infoMolt$M4.entry[2]/6); abline(v=infoMolt$M4.exit[2]/6)
plot(seq(0,60,(1/6)),gfpInten60); abline(v=infoMolt$hatch[3]/6); abline(v=infoMolt$M1.entry[3]/6); abline(v=infoMolt$M1.exit[3]/6); abline(v=infoMolt$M2.entry[3]/6); abline(v=infoMolt$M2.exit[3]/6); abline(v=infoMolt$M3.entry[3]/6); abline(v=infoMolt$M3.exit[3]/6); abline(v=infoMolt$M4.entry[3]/6); abline(v=infoMolt$M4.exit[3]/6)

#plot together
infoMoltc <- infoMolt[,5:13]-infoMolt$hatch+29 #set hatch t=29
gfpIntenc <- matrix(NA,nrow=nrow(infoMolt), ncol=402) #gfp for hatch at t=29
gfpIntenc[1,21:(length(gfpInten33)+20)] <- gfpInten33
gfpIntenc[2,8:(length(gfpInten34)+7)] <- gfpInten34
gfpIntenc[3,1:length(gfpInten60)] <- gfpInten60


gfpMoltc <- matrix(NA,nrow=nrow(infoMolt), ncol=402) #color hatch, M1, M2, M3 and M4 in blue
for (i in 1:nrow(infoMolt)){
  gfpMoltc[i,1:infoMoltc[i,1]] <- gfpIntenc[i,1:infoMoltc[i,1]] #hatch
  gfpMoltc[i,infoMoltc[i,2]:infoMoltc[i,3]] <- gfpIntenc[i,infoMoltc[i,2]:infoMoltc[i,3]] #M1
  gfpMoltc[i,infoMoltc[i,4]:infoMoltc[i,5]] <- gfpIntenc[i,infoMoltc[i,4]:infoMoltc[i,5]] #M2
  gfpMoltc[i,infoMoltc[i,6]:infoMoltc[i,7]] <- gfpIntenc[i,infoMoltc[i,6]:infoMoltc[i,7]] #M3
  gfpMoltc[i,infoMoltc[i,8]:infoMoltc[i,9]] <- gfpIntenc[i,infoMoltc[i,8]:infoMoltc[i,9]] #M4
}

gfpIntenc_mean <- na.omit(colMeans(gfpIntenc, na.rm = TRUE))


par(mfrow=c(2,1))
plot(seq(0,60,(1/6)),gfpInten34); abline(v=infoMolt$hatch[2]/6); abline(v=infoMolt$M1.entry[2]/6); abline(v=infoMolt$M1.exit[2]/6); abline(v=infoMolt$M2.entry[2]/6); abline(v=infoMolt$M2.exit[2]/6); abline(v=infoMolt$M3.entry[2]/6); abline(v=infoMolt$M3.exit[2]/6); abline(v=infoMolt$M4.entry[2]/6); abline(v=infoMolt$M4.exit[2]/6)
matplot(t(gfpIntenc), type='l', lty=1, col='grey', xaxt='n'); matplot(t(gfpMoltc), type='l', lty=2, col=alpha('blue',0.5), add=TRUE); lines(smooth.spline(gfpIntenc_mean, spar = 0.5), col='black'); axis(1, at = 1:403, labels = seq(0,67,1/6), cex.axis = 0.7)

