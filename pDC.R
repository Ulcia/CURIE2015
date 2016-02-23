setwd("/Users/urszulaczerwinska/Documents/CURIE2015/Data/DC_17.02.15")

pDC=read.table("HTSeq.txt", sep="\t", header=TRUE, row.names=1)
pDC=data.frame(pDC)
pDC.lg=log1p(pDC)
pDC.s=t(scale(t(pDC.lg), scale=FALSE))
summary(pDC.s)
#<<<<<<< HEAD
#=======
dim(pDC.s)
#>>>>>>> a8cfa4e

setwd("/Users/urszulaczerwinska/Documents/CURIE2015/Data/DC_17.02.15/pDC_ICA_all")
write.table(pDC.s, file= "pDCall_numerical.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
write.table(row.names(pDC.s), file= "pDCall_ids.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
write.table(colnames(pDC.s), file= "pDCall_samples.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
#<<<<<<< HEAD

#=======
write.table(pDC.s, file= "pDCall_full.txt", quote=FALSE, sep="\t", row.names=TRUE, col.names=TRUE)

pDC.s.nolog=t(scale(t(pDC), scale=FALSE))

pDC.an.nolog= merge(t(pDC.s.nolog), pDC_annot, by=0, all=TRUE) 
pDC.an <- merge(t(pDC.s), pDC_annot, by=0, all=TRUE) 
row.names(pDC.an)=pDC.an[,1]
row.names(pDC.an.nolog)=pDC.an.nolog[,1]

dim(pDC.an)
pDC.an=pDC.an[,c(2:14999, 15001, 15002)]
pDC.an.nolog=pDC.an.nolog[,c(2:14999, 15001, 15002)]

colnames(pDC.an)[14998:15000]=c("Tissue", "Subset", "Donor")
colnames(pDC.an.nolog)[14998:15000]=c("Tissue", "Subset", "Donor")

table(pDC.an$Subset)

PDC=subset(pDC.an, Subset=="PDC")[,1:14997]
MMAC=subset(pDC.an, Subset=="MMAC")[,1:14997]
BDCA1nDC=subset(pDC.an, Subset=="BDCA1nDC")[,1:14997]
CD14pDC=subset(pDC.an, Subset=="CD14pDC")[,1:14997]
BDCA1pDC=subset(pDC.an, Subset=="BDCA1pDC")[,1:14997]

PDC.s=t(scale(PDC,scale=FALSE))
MMAC.s=t(scale(MMAC,scale=FALSE))
BDCA1nDC.s=t(scale(BDCA1nDC,scale=FALSE))
CD14pDC.s=t(scale(CD14pDC,scale=FALSE))
BDCA1pDC.s=t(scale(BDCA1pDC,scale=FALSE))

setwd("/Users/urszulaczerwinska/Documents/CURIE2015/Data/DC_17.02.15/PDC")
write.table(PDC.s, file= "PDC.s_numerical.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
write.table(row.names(PDC.s), file= "PDC.s_ids.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
write.table(colnames(PDC.s), file= "PDC.s_samples.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
write.table(PDC.s, file= "PDC.s_full.txt", quote=FALSE, sep="\t", row.names=TRUE, col.names=TRUE)

setwd("/Users/urszulaczerwinska/Documents/CURIE2015/Data/DC_17.02.15/MMAC")
write.table(MMAC.s, file= "MMAC.s_numerical.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
write.table(row.names(MMAC.s), file= "MMAC.s_ids.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
write.table(colnames(MMAC.s), file= "MMAC.s_samples.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
write.table(MMAC.s, file= "MMAC.s_full.txt", quote=FALSE, sep="\t", row.names=TRUE, col.names=TRUE)

setwd("/Users/urszulaczerwinska/Documents/CURIE2015/Data/DC_17.02.15/BDCA1nDC")
write.table(BDCA1nDC.s, file= "BDCA1nDC.s_numerical.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
write.table(row.names(BDCA1nDC.s), file= "BDCA1nDC.s_ids.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
write.table(colnames(BDCA1nDC.s), file= "BDCA1nDC.s_samples.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
write.table(BDCA1nDC.s, file= "BDCA1nDC.s_full.txt", quote=FALSE, sep="\t", row.names=TRUE, col.names=TRUE)

setwd("/Users/urszulaczerwinska/Documents/CURIE2015/Data/DC_17.02.15/BDCA1pDC")
write.table(BDCA1pDC.s, file= "BDCA1pDC.s_numerical.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
write.table(row.names(BDCA1pDC.s), file= "BDCA1pDC.s_ids.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
write.table(colnames(BDCA1pDC.s), file= "BDCA1pDC.s_samples.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
write.table(BDCA1pDC.s, file= "BDCA1pDC.s_full.txt", quote=FALSE, sep="\t", row.names=TRUE, col.names=TRUE)

setwd("/Users/urszulaczerwinska/Documents/CURIE2015/Data/DC_17.02.15/CD14pDC")
write.table(CD14pDC.s, file= "CD14pDC.s_numerical.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
write.table(row.names(CD14pDC.s), file= "CD14pDC.s_ids.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
write.table(colnames(CD14pDC.s), file= "CD14pDC.s_samples.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
write.table(CD14pDC.s, file= "CD14pDC.s_full.txt", quote=FALSE, sep="\t", row.names=TRUE, col.names=TRUE)

summary(PDC.s)

#PCA
library(FactoMineR)
res.pca = PCA(pDC, scale.unit=TRUE, ncp=5, graph=T)
res.pca = PCA(pDC.an.nolog, scale.unit=TRUE, ncp=5, graph=F,quali.sup=c(14998,14999, 15000))
plot.PCA(res.pca, axes=c(1,2), choix="ind",habillage=14999)
#>>>>>>> a8cfa4e

setwd("/Users/urszulaczerwinska/Documents/CURIE2015/Data/DC_17.02.15")

pDC_annot=read.table("target_LABC.txt", sep="\t", header=TRUE, row.names=1)

summary(pDC_annot)
head(pDC_annot)

#<<<<<<< HEAD
#=======
setwd("/Users/urszulaczerwinska/Documents/CURIE2015/Data/DC_17.02.15/pDC_ICA_all")

#>>>>>>> a8cfa4e

pDC_2.S=read.table("S_pDCall_numerical.txt_2.num", sep="\t", header=FALSE, row.names=NULL)

row.names(pDC_2.S)=row.names(pDC.s)
colnames(pDC_2.S)=c("IC1","IC2")
pDC_2.S=pDC_2.S[,-3]
summary(pDC_2.S)
#<<<<<<< HEAD

pDC_2.S.s=t(scale(t(pDC_2.S),scale=FALSE))
#=======
dim(pDC_2.S)

hist(scale(pDC_2.S[,2]))

pDC_2.S.s=scale(pDC_2.S,scale=FALSE)
#>>>>>>> a8cfa4e
hist(pDC_2.S.s[,1])

hist(pDC_2.S.s[,2])

#<<<<<<< HEAD
#=======
pDC_2.S.s[,1][pDC_2.S.s[,1,drop=FALSE] < -2]
pDC_2.S.s[,1][pDC_2.S.s[,1] > 2]


pDC_2.S.s[,2][pDC_2.S.s[,2] < -2]
pDC_2.S.s[,2][pDC_2.S.s[,2] > 2]

pDC_IC1_2.S_m=names(pDC_2.S.s[,1][pDC_2.S.s[,1] < -3])
pDC_IC1_2.S_p=names(pDC_2.S.s[,1][pDC_2.S.s[,1] > 3])


pDC_IC2_2.S_m=names(pDC_2.S.s[,2][pDC_2.S.s[,2] < -2.5])
pDC_IC2_2.S_p=names(pDC_2.S.s[,2][pDC_2.S.s[,2] > 2.5])

write.table(c("pDC_IC1_2p", "na", t(pDC_IC1_2.S_p)), file="pDC_IC1_2p.gmt", quote=FALSE, sep="\t", eol = "\t",col.names=FALSE, row.names=FALSE)              
write.table(c("pDC_IC1_2m", "na", t(pDC_IC1_2.S_m)), file="pDC_IC1_2m.gmt", quote=FALSE, sep="\t", eol = "\t",col.names=FALSE, row.names=FALSE)              
write.table(c("pDC_IC2_2p", "na", t(pDC_IC2_2.S_p)), file="pDC_IC2_2p.gmt", quote=FALSE, sep="\t", eol = "\t",col.names=FALSE, row.names=FALSE)              
write.table(c("pDC_IC2_2m","na", t(pDC_IC2_2.S_m)), file="pDC_IC2_2m.gmt", quote=FALSE, sep="\t", eol = "\t",col.names=FALSE, row.names=FALSE)              

#########################
#>>>>>>> a8cfa4e
pDC_2.A=read.table("A_pDCall_numerical.txt_2.num", sep="\t", header=FALSE, row.names=NULL)
row.names(pDC_2.A)=colnames(pDC.s)
colnames(pDC_2.A)=c("IC1","IC2")
pDC_2.A=pDC_2.A[,-3]


summary(pDC_2.A)

#<<<<<<< HEAD
pDC_2.A.s=t(scale(t(pDC_2.A),scale=FALSE))
#=======
pDC_2.A.s=scale(pDC_2.A,scale=FALSE)
#>>>>>>> a8cfa4e
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

#HEAD
#=======
table(pDC_2.A.ann$Subset)


#>>>>>>> a8cfa4e
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


#<<<<<<< HEAD
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

#=======
x=pDC_2.A.ann[,1:2]
y=pDC_2.A.ann$Subset
featurePlot(x=x, y=y, plot="box",labels=c("Subset",""), main="Subset")

x=pDC_2.A.ann[,1:2]
y=pDC_2.A.ann$Tissue
featurePlot(x=x, y=y, plot="box", labels=c("Tissue",""), main="Tissue")

x=pDC_2.A.ann[,1:2]
y=pDC_2.A.ann$Donor
featurePlot(x=x, y=y, plot="box", labels=c("Donor",""), main="Donor")

####
x=pDC_2.A.ann[,1:2]
y=pDC_2.A.ann$Subset
featurePlot(x=x, y=y, plot="density",labels=c("Subset",""), main="Subset",auto.key=list(columns=5))

x=pDC_2.A.ann[,1:2]
y=pDC_2.A.ann$Tissue
featurePlot(x=x, y=y, plot="density", labels=c("Tissue",""), main="Tissue",auto.key=list(columns=2))

x=pDC_2.A.ann[,1:2]
y=pDC_2.A.ann$Donor
featurePlot(x=x, y=y, plot="density", labels=c("Donor",""), main="Donor",auto.key=list(columns=13))

#for continous varaibles
x=pDC_2.A.ann[,1:2]
y=pDC_2.A.ann$Subset
# scatterplot matrix
featurePlot(x=x, y=y, plot="ellipse",auto.key=list(columns=5))



############



dev.off()
pDC_6.S=read.table("S_pDCall_numerical.txt_6.num", sep="\t", header=FALSE, row.names=NULL)

row.names(pDC_6.S)=row.names(pDC.s)
colnames(pDC_6.S)=c("IC1","IC2", "IC3", "IC4", "IC5", "IC6")
pDC_6.S=pDC_6.S[,-7]
summary(pDC_6.S)
dim(pDC_6.S)


pDC_6.S.s=scale(pDC_6.S,scale=FALSE)
par(mfrow=c(2,3))
hist(pDC_6.S.s[,1])

hist(pDC_6.S.s[,2])
hist(pDC_6.S.s[,3])
hist(pDC_6.S.s[,4])
hist(pDC_6.S.s[,5])
hist(pDC_6.S.s[,6])

pDC_6.S.s[,1][pDC_6.S.s[,1] < -2]
pDC_6.S.s[,1][pDC_6.S.s[,1] > 2]

pDC_6.S.s[,2][pDC_6.S.s[,2] < -2]
pDC_6.S.s[,2][pDC_6.S.s[,2] > 2]

pDC_6.S.s[,3][pDC_6.S.s[,3] < -2]
pDC_6.S.s[,3][pDC_6.S.s[,3] > 2]

pDC_6.S.s[,4][pDC_6.S.s[,4] < -2]
pDC_6.S.s[,4][pDC_6.S.s[,4] > 2]

pDC_6.S.s[,5][pDC_6.S.s[,5] < -2]
pDC_6.S.s[,5][pDC_6.S.s[,5] > 2]

pDC_6.S.s[,6][pDC_6.S.s[,6] < -2]
pDC_6.S.s[,6][pDC_6.S.s[,6] > 2]


write.table(pDC_6.S.s[,1][pDC_6.S.s[,1] < -2], file="pDC_IC1_6m.rnk", quote=FALSE, sep="\t", col.names=FALSE, row.names=TRUE)              

write.table(pDC_6.S.s[,1][pDC_6.S.s[,1] > 2], file="pDC_IC1_6p.rnk", quote=FALSE, sep="\t", col.names=FALSE, row.names=TRUE)              
write.table(pDC_6.S.s[,2][pDC_6.S.s[,2] < -2], file="pDC_IC2_6m.rnk", quote=FALSE, sep="\t", col.names=FALSE, row.names=TRUE)              
write.table(pDC_6.S.s[,2][pDC_6.S.s[,2] > 2], file="pDC_IC2_6p.rnk", quote=FALSE, sep="\t", col.names=FALSE, row.names=TRUE)              
write.table(pDC_6.S.s[,3][pDC_6.S.s[,3] < -2], file="pDC_IC3_6m.rnk", quote=FALSE, sep="\t", col.names=FALSE, row.names=TRUE)              
write.table(pDC_6.S.s[,3][pDC_6.S.s[,3] > 2], file="pDC_IC3_6p.rnk", quote=FALSE, sep="\t", col.names=FALSE, row.names=TRUE)              
write.table(pDC_6.S.s[,4][pDC_6.S.s[,4] < -2], file="pDC_IC4_6m.rnk", quote=FALSE, sep="\t", col.names=FALSE, row.names=TRUE)              
write.table(pDC_6.S.s[,4][pDC_6.S.s[,4] > 2], file="pDC_IC4_6p.rnk", quote=FALSE, sep="\t", col.names=FALSE, row.names=TRUE)              

write.table(pDC_6.S.s[,5][pDC_6.S.s[,5] < -2], file="pDC_IC5_6m.rnk", quote=FALSE, sep="\t", col.names=FALSE, row.names=TRUE)              
write.table(pDC_6.S.s[,5][pDC_6.S.s[,5] > 2], file="pDC_IC5_6p.rnk", quote=FALSE, sep="\t", col.names=FALSE, row.names=TRUE)              
write.table(pDC_6.S.s[,6][pDC_6.S.s[,6] < -2], file="pDC_IC6_6m.rnk", quote=FALSE, sep="\t", col.names=FALSE, row.names=TRUE)              
write.table(pDC_6.S.s[,6][pDC_6.S.s[,6] > 2], file="pDC_IC6_6p.rnk", quote=FALSE, sep="\t", col.names=FALSE, row.names=TRUE)              

#