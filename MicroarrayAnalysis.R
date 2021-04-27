library(limma)
library(affy)

setwd("E:/MicroArrayData/DBT core proposal/Microarray Processing/")

targets=read.csv("Targets.csv")
raw = read.maimages(targets, source="agilent",green.only=TRUE)
raw_BGcorrected = backgroundCorrect(raw, method="normexp", offset=16)
raw_BGandNormalized = normalizeBetweenArrays(raw_BGcorrected,method="quantile")
raw_aver = avereps(raw_BGandNormalized,ID=raw_BGandNormalized$genes$ProbeName)

boxplot(log(as.matrix(raw_BGcorrected)),las=2,ylab="Log2(Intensity)")

boxplot(as.matrix(raw_BGandNormalized), las=2, ylab = "Log2(Intensity)")

f = factor(targets$Category)
design = model.matrix(~f)

design = cbind(R = c(1,1,1,0,0,0,0,0),   
               S = c(0,0,0,1,1,1,0,0),    
               CR = c(0,0,0,0,0,0,1,0),    
               CS = c(0,0,0,0,0,0,0,1)) 

fit = lmFit(raw_aver, design)
contrastMatrix = makeContrasts("CR-CS", levels=design)

fit2 = contrasts.fit(fit, contrastMatrix)
fit2 = eBayes(fit2)

topTable(fit2, coef = "CR-CS")


sig = length(which(topTable(fit2, coef = "CR-CS",number=15744)[,15]<0.05))

signif = topTable(fit2, coef = "CR-CS",number=sig)
upregulated = signif[which(signif[,11]>0),]
downregulated = signif[which(signif[,11]<0),]

#Normalized Gene Expression for GIMME
write.csv(as.matrix(raw_BGandNormalized),"Normalized.csv")

#UpregulatedandDownregulated
write.csv(upregulated, "CRCS_Upre.csv")
write.csv(downregulated, "CRCS_Downre.csv")

