setwd("/Users/urszulaczerwinska/Documents/CURIE2015/Data/DC_17.02.15")

pDC=read.table("HTSeq.txt", sep="\t", header=TRUE, row.names=1)
pDC=data.frame(pDC)
pDC.lg=log1p(pDC)
pDC.s=t(scale(t(pDC.lg), scale=FALSE))
summary(pDC.s)

setwd("/Users/urszulaczerwinska/Documents/CURIE2015/Data/DC_17.02.15/pDC_ICA_all")
write.table(pDC.s, file= "pDCall_numerical.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
write.table(row.names(pDC.s), file= "pDCall_ids.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
write.table(colnames(pDC.s), file= "pDCall_samples.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)


setwd("/Users/urszulaczerwinska/Documents/CURIE2015/Data/DC_17.02.15")

pDC_annot=read.table("target_LABC.txt", sep="\t", header=TRUE, row.names=1)

summary(pDC_annot)
head(pDC_annot)


pDC_2.S=read.table("S_pDCall_numerical.txt_2.num", sep="\t", header=FALSE, row.names=NULL)

row.names(pDC_2.S)=row.names(pDC.s)
colnames(pDC_2.S)=c("IC1","IC2")
pDC_2.S=pDC_2.S[,-3]
summary(pDC_2.S)

pDC_2.S.s=t(scale(t(pDC_2.S),scale=FALSE))
hist(pDC_2.S.s[,1])

hist(pDC_2.S.s[,2])

pDC_2.A=read.table("A_pDCall_numerical.txt_2.num", sep="\t", header=FALSE, row.names=NULL)
row.names(pDC_2.A)=colnames(pDC.s)
colnames(pDC_2.A)=c("IC1","IC2")
pDC_2.A=pDC_2.A[,-3]


summary(pDC_2.A)

pDC_2.A.s=t(scale(t(pDC_2.A),scale=FALSE))
hist(pDC_2.A.s[,1])

hist(pDC_2.A.s[,2])

pDC_2.A.s[,1][pDC_2.A.s[,1] < -0.1]
pDC_2.A.s[,1][pDC_2.A.s[,1] > 0.1]


pDC_2.A.s[,2][pDC_2.A.s[,2] < -0.1]
pDC_2.A.s[,2][pDC_2.A.s[,2] > 0.1]

head(pDC_2.A.s)


pDC_2.A.ann <- merge(pDC_2.A.s, pDC_annot, by=0, all=TRUE) 
summary(pDC_2.A.ann)
row.names(pDC_2.A.ann)=pDC_2.A.ann[,1]
head(pDC_2.A.ann)
dim(pDC_2.A.ann)

pDC_2.A.ann=pDC_2.A.ann[,2:12]

library(caret)
library(AppliedPredictiveModeling)
transparentTheme(trans = .4)

theme1 <- trellis.par.get()
theme1$plot.symbol$col = rgb(.2, .2, .2, .4)
theme1$plot.symbol$pch = 16
theme1$plot.line$col = rgb(1, 0, 0, .7)
theme1$plot.line$lwd <- 2
trellis.par.set(theme1)


x=pDC_2.A.ann$IC1
y=pDC_2.A.ann$Subset
featurePlot(x=x, y=y, plot="box")

x=pDC_2.A.ann$IC2
y=pDC_2.A.ann$Subset
featurePlot(x=x, y=y, plot="box")

x=pDC_2.A.ann$IC2
y=pDC_2.A.ann$Tissue
featurePlot(x=x, y=y, plot="box")

x=pDC_2.A.ann$IC1
y=pDC_2.A.ann$Tissue
featurePlot(x=x, y=y, plot="box")

