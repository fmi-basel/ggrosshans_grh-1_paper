writeLines(capture.output(sessionInfo()), "Heatmap_L1phen_sessionInfo.txt")

library('R.matlab')
library(NMF)
setwd("/tungstenfs/scratch/ggrossha/meeumilo/scripts/ggrosshans_grh-1_paper")

file <- readMat('Data/2018-03-27 grh-1_highAuxRange.mat')
valid <- as.factor(file$valid)
annotation <- unlist(file$samples)
annotation <- ordered(annotation, levels=c("ctr", "60.9 nM", "243 nM", "975 nM", "3.9 uM", "15.6 uM", '62.5 uM', '250 uM', '1 mM')) #order the annotation to get right order in boxplot

##- duration stages
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
df_250 <- df[which(df$condition=='250 uM' & df$valid==1),] 

## plot L1 phenotype
L1phen_df_aux <- df_250[order(df_250$M1_in),]
L1phen <- as.numeric(rownames(L1phen_df_aux))[1:20]
myCol <- colorRampPalette(c('black','white'))(100)
aheatmap(t(file$Xc[,L1phen]), Rowv = NA, Colv = NA, breaks=seq(-1,4,length.out = 51), color = myCol, main='L1 phenotype')
