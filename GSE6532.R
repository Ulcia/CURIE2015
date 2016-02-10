setwd("/Users/urszulaczerwinska/Documents/CURIE2015/Data/Biton.Zinovyev/GSE6532/ICAtest/8ics/")
IC1<- read.table("./input_GSEA1.rnk", sep="\t", header=F)

IC2<- read.table("./input_GSEA2.rnk", sep="\t", header=F)
IC3<- read.table("./input_GSEA3.rnk", sep="\t", header=F)
IC4<- read.table("./input_GSEA4.rnk", sep="\t", header=F)
IC5<- read.table("./input_GSEA5.rnk", sep="\t", header=F)
IC6<- read.table("./input_GSEA6.rnk", sep="\t", header=F)
IC7<- read.table("./input_GSEA7.rnk", sep="\t", header=F)
IC8<- read.table("./input_GSEA8.rnk", sep="\t", header=F)
IC8Bitton=read.table("/Users/urszulaczerwinska/Documents/CURIE2015/Data/Biton.Zinovyev/IC8_immune_projections",row.names = 1, sep="\t", header=T)
intersect <- function(x, y) y[match(x, y, nomatch = 0)]

summary(IC8Bitton)
head(IC8Bitton)
dev.off()
par(mfrow=c(2,8))

boxplot(IC3[,2])



hist(IC1[,2], col="black")
barplot(IC1[,2], main="barplot IC1")


hist(IC2[,2], col="black")
barplot(IC2[,2], main="barplot IC2")

hist(IC3[,2], col="black")
barplot(IC3[,2], main="barplot IC3")


hist(IC4[,2], col="black")
barplot(IC4[,2], main="barplot IC4")

hist(IC5[,2], col="black")
barplot(IC5[,2], main="barplot IC5")

hist(IC6[,2], col="black")
barplot(IC6[,2], main="barplot IC6")

hist(IC7[,2], col="black")
barplot(IC7[,2], main="barplot IC7")

hist(IC8[,2], col="black")
barplot(IC8[,2], main="barplot IC8")


myICcorr(IC1)
myICcorr(IC2)
#myICcorr(IC3,scale=TRUE)
myICcorr(IC3)

myICcorr(IC4)
myICcorr(IC5)
myICcorr(IC6)
myICcorr(IC7)
myICcorr(IC8)

B=cbind(myICcorr(IC1),
myICcorr(IC2),
#myICcorr(IC3,scale=TRUE)
myICcorr(IC3),

myICcorr(IC4),
myICcorr(IC5),
myICcorr(IC6),
myICcorr(IC7),
myICcorr(IC8))

corr.mat=cbind(IC8Bitton.i[b,1],B)
test.matrix.cor=rcorr(corr.mat)
test.matrix.cor$P[1,1]=0
test.matrix.cor$P[2,2]=0
test.matrix.cor$P[3,3]=0
test.matrix.cor$P[4,4]=0
test.matrix.cor$P[5,5]=0
test.matrix.cor$P[6,6]=0
test.matrix.cor$P[7,7]=0
test.matrix.cor$P[8,8]=0


test.matrix.cor.f=flattenCorrMatrix(test.matrix.cor$r, test.matrix.cor$P)
corrplot(test.matrix.cor$r[,1], type="lower", order="original", tl.col="black", tl.srt=45, col=colmy(300))

data.frame(cbind(test.matrix.cor$r[2:9,1],test.matrix.cor$P[2:9,1]))

dataset=IC1

myICcorr= function(dataset,method="",scale=FALSE) {
  
datasetu <- subset(dataset, !duplicated(dataset[,1]))
#IC2u <- subset(IC2, !duplicated(IC2[,1]))
#IC3u <- subset(IC3, !duplicated(IC3[,1]))
#IC6u <- subset(IC6, !duplicated(IC6[,1]))

#write.table(IC8u, file= "input_GSEA8_u.rnk", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)



row.names(datasetu)=datasetu[,1]
#row.names(IC2u)=IC2u[,1]
#row.names(IC8u)=IC8u[,1]

myrowsImm=intersect(row.names(datasetu),row.names(IC8Bitton))
#myrowsImm=intersect(row.names(IC1u),row.names(IC4Bitton))
#myrowsImm=intersect(row.names(IC8u),row.names(IC4Bitton))



datasetu.i=datasetu[myrowsImm,,drop=FALSE]
IC8Bitton.i=IC8Bitton[myrowsImm,,drop=FALSE]

#IC2u.i=IC2u[myrowsImm,,drop=FALSE]
#IC4Bitton.i=IC4Bitton[myrowsImm,,drop=FALSE]

#IC8u.i=IC8u[myrowsImm,,drop=FALSE]
#IC4Bitton.i=IC4Bitton[myrowsImm,,drop=FALSE]



#row.names(IC4Bitton.i)
#IC2u[myrowsImm,2]
#IC4Bitton[myrowsImm,1]
o=order(row.names(datasetu.i))
b=order(row.names(IC8Bitton.i))

#head(IC8u.i[o,])
#head(IC4Bitton.i[o,, drop=FALSE])


if (scale==TRUE)
{
if (method == "spearman"){
cor.test(scale(datasetu.i[o,2]),scale(IC8Bitton.i[b,1]), method="spearman")
}
else
  {cor.test(scale(datasetu.i[o,2]),scale(IC8Bitton.i[b,1]))}

}
else{
  if (method == "spearman"){
    cor.test(datasetu.i[o,2],IC8Bitton.i[b,1], method="spearman")
  }
  else
  {cor.test(datasetu.i[o,2],IC8Bitton.i[b,1])}
  
}

#return (tocorr=datasetu.i[o,2])
}


row.names(IC6u)=IC6u[,1]

head(IC6u)

dataset=IC3
datasetu <- subset(dataset, !duplicated(dataset[,1]))
datasetu =data.frame(datasetu)
datasetu=na.omit(datasetu)
datasetu.o=datasetu[order(datasetu[,2]),]
write.table(datasetu.o, file= "input_GSEA3_u.rnk", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

sd(scale(datasetu[,2]))
3*sd(datasetu[,2])
mean(datasetu[,2])
mean(datasetu[,2])+3*sd(datasetu[,2])
mean(datasetu[,2])-3*sd(datasetu[,2])

datasetu.f.plus=datasetu[datasetu$V2 > 3.339,]

datasetu.f.plus=datasetu.o[scale(datasetu.o$V2) > 2.5,]
datasetu.f.min=datasetu.o[scale(datasetu.o$V2)< -2.5,]
datasetu.f=datasetu.o[abs(scale(datasetu.o$V2))> 2.5,]
dim(datasetu.f.min)
write.table(datasetu.f.plus, file= "input_GSEA3_plus.rnk", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
write.table(datasetu.f.min, file= "input_GSEA3_min.rnk", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
write.table(datasetu.f, file= "input_GSEA3_filtered.rnk", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

setwd("/Users/urszulaczerwinska/Documents/CURIE2015/Data/Biton.Zinovyev/GSE6532/ICAtest/8ics/IC3_BreastBCR_min")


IC1.IC3<- read.table("./ranked_IC1_IC3_breastBCR_min.rnk", sep="\t", header=T, row.names = 1)

IC2.IC3<- read.table("./ranked_IC2_IC3_breastBCR_min.rnk", sep="\t", header=T, row.names = 1)
IC3.IC3<- read.table("./ranked_IC3_IC3_breastBCR_min.rnk", sep="\t", header=T, row.names = 1)
IC4.IC3<- read.table("./ranked_IC4_IC3_breastBCR_min.rnk", sep="\t", header=T, row.names = 1)
IC5.IC3<- read.table("./ranked_IC5_IC3_breastBCR_min.rnk", sep="\t", header=T, row.names = 1)


IC1.IC3.s=scale(IC1.IC3)
IC2.IC3.s=scale(IC2.IC3)
IC3.IC3.s=scale(IC3.IC3)
IC4.IC3.s=scale(IC4.IC3)
IC5.IC3.s=scale(IC5.IC3)


IC1.n<- IC1.IC3.s[ order(row.names(IC1.IC3.s)), ]
IC2.n<- IC2.IC3.s[ order(row.names(IC2.IC3.s)), ]
cor.test(IC1.n[sig.rows],IC2.n[sig.rows])




IC1.IC3_min=IC1.IC3.s[IC1.IC3.s< -1,, drop=FALSE]
IC2.IC3_min=IC2.IC3.s[IC2.IC3.s> 2,, drop=FALSE]
IC2.IC3_plus=IC2.IC3.s[IC2.IC3.s >2,, drop=FALSE]
IC3.IC3_plus=IC3.IC3.s[IC3.IC3.s < -2,, drop=FALSE]
IC3.IC3_plus=IC3.IC3.s[IC3.IC3.s >2,, drop=FALSE]
IC4.IC3_plus=IC4.IC3.s[IC4.IC3.s < -2,, drop=FALSE]
IC5.IC3_plus=IC5.IC3.s[IC5.IC3.s < -1,, drop=FALSE]

IC1.IC3_min=IC1.IC3.s[IC1.IC3.s> 1,, drop=FALSE]


length(which(row.names(IC1.IC3.s)%in%names(m$B)==TRUE))#2
length(which(row.names(IC1.IC3.s)%in%names(m$T)==TRUE))#2
length(which(row.names(IC1.IC3.s)%in%names(m$LYMPHS)==TRUE)) #2
length(which(row.names(IC1.IC3.s)%in%names(m$GRANS)==TRUE)) #2
length(which(row.names(IC1.IC3.s)%in%names(m$CD8)==TRUE)) #2



length(which(row.names(IC5.IC3_plus)%in%names(m$B)==TRUE))#2
length(which(row.names(IC5.IC3_plus)%in%names(m$T)==TRUE))#2
length(which(row.names(IC5.IC3_plus)%in%names(m$LYMPHS)==TRUE)) #2
length(which(row.names(IC5.IC3_plus)%in%names(m$GRANS)==TRUE)) #2
length(which(row.names(IC5.IC3_plus)%in%names(m$CD8)==TRUE)) #2


length(which(row.names(IC1.IC3_min)%in%names(m$B)==TRUE))#2



write.table(IC1.IC3_min, file= "IC1.IC3_min_wv.rnk", quote=FALSE, sep="\t", row.names=TRUE, col.names=FALSE)
write.table(IC2.IC3_min, file= "IC2.IC3_min._wv.rnk", quote=FALSE, sep="\t", row.names=TRUE, col.names=FALSE)
write.table(IC2.IC3_plus, file= "IC2.IC3_plus_wv.rnk", quote=FALSE, sep="\t", row.names=TRUE, col.names=FALSE)
write.table(IC3.IC3_plus, file= "IC3.IC3_plus_wv.rnk", quote=FALSE, sep="\t", row.names=TRUE, col.names=FALSE)
write.table(IC4.IC3_plus, file= "IC4.IC3_plus_wv.rnk", quote=FALSE, sep="\t", row.names=TRUE, col.names=FALSE)
write.table(IC5.IC3_plus, file= "IC5.IC3_plus_wv.rnk", quote=FALSE, sep="\t", row.names=TRUE, col.names=FALSE)



hist(IC1.IC3, col="black")
barplot(IC2.IC3, main="barplot IC1")

library(CellMix)
m <- MarkerList("Palmer")
summary(m)
c("PALMER Bcell",names(m$B))

setwd("/Users/urszulaczerwinska/Documents/CURIE2015/Data/Biton.Zinovyev/GSE6532/ICAtest/8ics/IC3_BreastBCR_min/Immune_sig_gmt")

write.table(t(c("PALMER Bcell", "na",names(m$B))), file= "PALMER_Bcell.gmt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

exprs(m)
write.table(t(c("PALMER CD8", "na",names(m$CD8))), file= "PALMER_CD8.gmt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

names(m$GRANS)
write.table(t(c("PALMER GRANS", "na",names(m$GRANS))), file= "PALMER_GRANS.gmt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

names(m$LYMPHS)

write.table(t(c("PALMER LYMPH", "na",names(m$LYMPHS))), file= "PALMER_LYMPHS.gmt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

names(m$T)
write.table(t(c("PALMER Tcell", "na",names(m$T))), file= "PALMER_T.gmt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

IRIS <- MarkerList("IRIS")
summary(IRIS)
IRIS.h <- convertIDs(IRIS, 'SYMBOL', verbose=2)
summary(IRIS.h)



IC1.IC3_plus=IC1.IC3.s[IC1.IC3.s< -2,, drop=FALSE]
IC2.IC3_plus=IC2.IC3.s[IC2.IC3.s < -1,, drop=FALSE]
IC3.IC3_plus=IC3.IC3.s[IC3.IC3.s > 1,, drop=FALSE]
IC3.IC3_plus=IC3.IC3.s[IC3.IC3.s < -2,, drop=FALSE]
IC4.IC3_plus=IC4.IC3.s[IC4.IC3.s > 2,, drop=FALSE]
IC5.IC3_plus=IC5.IC3.s[IC5.IC3.s > 2,, drop=FALSE]


length(which(row.names(IC1.IC3.s)%in%names(IRIS.h$B)==TRUE))#2
length(which(row.names(IC1.IC3.s)%in%names(IRIS.h$T)==TRUE))#2
length(which(row.names(IC1.IC3.s)%in%names(IRIS.h$Dendritic)==TRUE)) #2
length(which(row.names(IC1.IC3.s)%in%names(IRIS.h$NK)==TRUE)) #2
length(which(row.names(IC1.IC3.s)%in%names(IRIS.h$Monocyte)==TRUE)) #2
length(which(row.names(IC1.IC3.s)%in%names(IRIS.h$Neutrophil)==TRUE)) #2



length(which(row.names(IC2.IC3_plus)%in%names(IRIS.h$B)==TRUE))#2
length(which(row.names(IC2.IC3_plus)%in%names(IRIS.h$T)==TRUE))#2
length(which(row.names(IC2.IC3_plus)%in%names(IRIS.h$Dendritic)==TRUE)) #2
length(which(row.names(IC2.IC3_plus)%in%names(IRIS.h$NK)==TRUE)) #2
length(which(row.names(IC2.IC3_plus)%in%names(IRIS.h$Monocyte)==TRUE)) #2
length(which(row.names(IC2.IC3_plus)%in%names(IRIS.h$Neutrophil)==TRUE)) #2


setwd("/Users/urszulaczerwinska/Documents/CURIE2015/Data/Biton.Zinovyev/GSE6532/ICAtest/8ics/IC3_BreastBCR_min/Immune_sig_gmt")

write.table(t(c("IRIS Bcell", "na",names(IRIS.h$B))), file= "IRIS_Bcell.gmt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

names(IRIS.h$T)
write.table(t(c("IRIS Tcell", "na",names(IRIS.h$T))), file= "IRIS_Tcell.gmt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

names(IRIS.h$Dendritic)
write.table(t(c("IRIS Dendritic", "na",names(IRIS.h$Dendritic))), file= "IRIS_Dendritic.gmt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

names(IRIS.h$NK)

write.table(t(c("IRIS NK", "na",names(IRIS.h$NK))), file= "IRIS_NK.gmt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

names(IRIS.h$Monocyte)
write.table(t(c("IRIS Monocyte", "na",names(IRIS.h$Monocyte))), file= "IRIS_Monocyte.gmt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

names(IRIS.h$Neutrophil)
write.table(t(c("IRIS Neutrophil", "na",names(IRIS.h$Neutrophil))), file= "IRIS_Neutrophil.gmt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

names(IRIS.h$Lymphoid)
write.table(t(c("IRIS Lymphoid", "na",names(IRIS.h$Lymphoid))), file= "IRIS_Lymphoid.gmt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

names(IRIS.h$Myeloid)
write.table(t(c("IRIS Myeloid", "na",names(IRIS.h$Myeloid))), file= "IRIS_Myeloid.gmt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

names(IRIS.h$Multiple)
write.table(t(c("IRIS Multiple", "na",names(IRIS.h$Multiple))), file= "IRIS_Multiple.gmt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)


Ha <- MarkerList("HaemAtlas")
summary(Ha)
Ha.h <- convertIDs(Ha, 'SYMBOL', verbose=2)
summary(Ha.h)



IC1.IC3_min=IC1.IC3.s[IC1.IC3.s< -1,, drop=FALSE]
IC2.IC3_min=IC2.IC3.s[IC2.IC3.s< -2,, drop=FALSE]
IC2.IC3_plus=IC2.IC3.s[IC2.IC3.s >2,, drop=FALSE]
IC3.IC3_plus=IC3.IC3.s[IC3.IC3.s > 1,, drop=FALSE]
IC3.IC3_plus=IC3.IC3.s[IC3.IC3.s >2,, drop=FALSE]
IC4.IC3_plus=IC4.IC3.s[IC4.IC3.s > 1,, drop=FALSE]
IC5.IC3_plus=IC5.IC3.s[IC5.IC3.s < -2,, drop=FALSE]


length(which(row.names(IC1.IC3.s)%in%Ha.h$'B-CD19'==TRUE))#2
length(which(row.names(IC1.IC3.s)%in%Ha.h$Erythroblast==TRUE))#2
length(which(row.names(IC1.IC3.s)%in%Ha.h$'Granulocyte-CD66b'==TRUE)) #2
length(which(row.names(IC1.IC3.s)%in%Ha.h$Megakaryocyte==TRUE)) #2
length(which(row.names(IC1.IC3.s)%in%Ha.h$'Monocyte-CD14'==TRUE)) #2
length(which(row.names(IC1.IC3.s)%in%Ha.h$'NK-CD56'==TRUE)) #2


write.table(t(c("Haem Atlas B-CD19", "na",Ha.h$'B-CD19')), file= "H_atlas_B-CD19.gmt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

Ha.h$Erythroblast
write.table(t(c("Haem Atlas Erythroblast", "na",Ha.h$Erythroblast)), file= "H_atlas_Erythroblast.gmt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

Ha.h$'Granulocyte-CD66b'
write.table(t(c("Haem Atlas Granulocyte-CD66b", "na",Ha.h$'Granulocyte-CD66b')), file= "H_atlas_Granulocyte-CD66b.gmt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

Ha.h$Megakaryocyte

write.table(t(c("Haem Atlas Megakaryocyte", "na",Ha.h$Megakaryocyte)), file= "H_atals_Megakaryocyte.gmt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

Ha.h$'Monocyte-CD14'
write.table(t(c("Haem Atlas Monocyte-CD14", "na",Ha.h$'Monocyte-CD14')), file= "H_atlas_Monocyte-CD14.gmt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

Ha.h$'Monocyte-CD14'
write.table(t(c("Haem Atlas Monocyte-CD14", "na",Ha.h$'Monocyte-CD14')), file= "H_atlas_Monocyte-CD14.gmt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

Ha.h$'NK-CD56'
write.table(t(c("Haem Atlas NK-CD56", "na",Ha.h$'NK-CD56')), file= "H_atlas_NK-CD56.gmt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)



G <- MarkerList("Grigoryev")
summary(G)
G.h <- convertIDs(G, 'SYMBOL', verbose=2)
summary(G.h)

G.h$'B cells'
write.table(t(c("Grigoryev B cells", "na",G.h$'B cells')), file= "Grigoryev_Bcells .gmt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

G.h$Eosinophils
write.table(t(c("Grigoryev Eosinophils", "na",G.h$Eosinophils)), file= "Grigoryev_Eosinophils.gmt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

G.h$'Monocytes'   
write.table(t(c("Grigoryev Monocyte", "na",G.h$'Monocytes')), file= "Grigoryev_Monocyte.gmt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

G.h$Neutrophils  
write.table(t(c("Grigoryev Neutrophils", "na",G.h$Neutrophils)), file= "Grigoryev_Neutrophils.gmt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

G.h$'NK'  
write.table(t(c("Grigoryev NK", "na",G.h$'NK')), file= "Grigoryev_NK.gmt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

G.h$'T cells'  
write.table(t(c("Grigoryev T cells", "na",G.h$'T cells')), file= "Grigoryev_Tcells.gmt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)


A<- MarkerList("Abbas")
summary(A)
A.h <- convertIDs(A, 'SYMBOL', verbose=2)
summary(A.h)



names(A.h$'Th act')
write.table(t(c("Abbas Th", "na",names(A.h$'Th act'))), file= "Abbas_Th.gmt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

names(A.h$'Tc act')
write.table(t(c("Abbas Tc", "na",c(names(A.h$'Tc act'),names(A.h$'Tc')))), file= "Abbas_Tc.gmt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

names(A.h$'B act')
write.table(t(c("Abbas B", "na",names(A.h$'B act'))), file= "Abbas_B.gmt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

names(A.h$'PC')
write.table(t(c("Abbas PC", "na",names(A.h$'PC'))), file= "Abbas_PC.gmt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

write.table(t(c("Abbas NK", "na",c(names(A.h$'NK act'),names(A.h$'NK')))), file= "Abbas_NK.gmt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

write.table(t(c("Abbas mono", "na",c(names(A.h$'mono act'),names(A.h$'mono')))), file= "Abbas_mono.gmt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

write.table(t(c("Abbas DC", "na",c(names(A.h$'DC act'),names(A.h$'DC')))), file= "Abbas_DC.gmt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
write.table(t(c("Abbas neutro", "na",names(A.h$'neutro'))), file= "Abbas_neutro.gmt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)


write.table(t(c("Bcell_CellMix", "na",unique(c(names(IRIS.h$B),names(m$B),Ha.h$'B-CD19',G.h$'B cells',names(A.h$'B act'))))), file= "Bcell.gmt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
write.table(t(c("Tcell_CellMix", "na",unique(c(names(IRIS.h$T),names(m$T),G.h$'T cells',names(A.h$'Tc act'),names(A.h$'Th act'),names(A.h$'Tc'))))), file= "Tcell.gmt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
write.table(t(c("mono_CellMix", "na",unique(c(names(IRIS.h$Monocyte),Ha.h$'Monocyte-CD14',G.h$'Monocytes',names(A.h$'mono act'))))), file= "mono.gmt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
write.table(t(c("NK_CellMix", "na",unique(c(names(IRIS.h$NK),Ha.h$'NK-CD56',G.h$'NK',names(A.h$'NK act'),names(A.h$'NK'))))), file= "NK.gmt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

#--------GSEA
cosmos=format_from_gmt("/Users/urszulaczerwinska/Documents/CURIE2015/ICA_Arnau/ICA_script/msigdb.v4.0.symbols_KEGG.gmt")
c7=format_from_gmt_2("/Users/urszulaczerwinska/Documents/CURIE2015/GSEA/msigdb_v5.0_files_to_download_locally/msigdb_v5.0_GMTs/c7.all.v5.0.symbols.gmt")
file="./Abbas.gmt"
file2="./Tcell.gmt"
setwd("/Users/urszulaczerwinska/Documents/CURIE2015/Data/Biton.Zinovyev/GSE6532/ICAtest/8ics/IC3_BreastBCR_min/Immune_sig_gmt/")
Abbas<-format_from_gmt_2("./Abbas.gmt")
Tcell<-format_from_gmt_2(file2)
Bcell<-format_from_gmt_2("./Bcell.gmt")
IRIS<-format_from_gmt_2("./IRIS.gmt")
mono<-format_from_gmt_2("./mono.gmt")
NK<-format_from_gmt_2("./NK.gmt")
Grigo<-format_from_gmt_2("./Grigoryev.gmt")
Ha<-format_from_gmt_2("./Haem_atlas.gmt")
Palmer<-format_from_gmt_2("./Palmer.gmt")

LM22_gmt=format_from_gmt_2("./LM22.gmt")
Freeman=format_from_gmt_2("./Freeman.gmt")

genes_test

Example<-enrichment(genes_test, min_module_size = 10, threshold = 0.05, maps = list(Abbas=Abbas, Tcell=Tcell))


Example[[1]]
#------------


IC1.IC3_min=IC1.IC3.s[IC1.IC3.s< -1,, drop=FALSE]
IC2.IC3_plus=IC2.IC3.s[IC2.IC3.s > 1,, drop=FALSE]
IC3.IC3_plus=IC3.IC3.s[IC3.IC3.s > 1,, drop=FALSE]
IC4.IC3_plus=IC4.IC3.s[IC4.IC3.s > 1,, drop=FALSE]
IC5.IC3_plus=IC5.IC3.s[IC5.IC3.s > 1,, drop=FALSE]
setwd("/Users/urszulaczerwinska/Documents/CURIE2015/Data/Biton.Zinovyev/GSE6532/Cytoscape/")
write.table(row.names(IC1.IC3_min),file="IC1.IC3_tr1.txt",sep="\t", col.names = FALSE,row.names=FALSE, quote =FALSE)
write.table(row.names(IC2.IC3_plus),file="IC2.IC3_tr1.txt",sep="\t", col.names = FALSE,row.names=FALSE, quote =FALSE)
write.table(row.names(IC3.IC3_plus),file="IC3.IC3_tr1.txt",sep="\t", col.names = FALSE,row.names=FALSE, quote =FALSE)
write.table(row.names(IC4.IC3_plus),file="IC4.IC3_tr1.txt",sep="\t", col.names = FALSE,row.names=FALSE, quote =FALSE)
write.table(row.names(IC5.IC3_plus),file="IC5.IC3_tr1.txt",sep="\t", col.names = FALSE,row.names=FALSE, quote =FALSE)

All_sign=c(row.names(IC1.IC3_min),row.names(IC2.IC3_plus),row.names(IC3.IC3_plus),row.names(IC4.IC3_plus),row.names(IC5.IC3_plus))
All_sign.cl=subset(All_sign, !duplicated(All_sign))
All_sign.cl.m=data.tam.t.n.u[All_sign.cl,]
dim(All_sign.cl.m)
All_sign.cl.m.c=t(scale(t(All_sign.cl.m),scale = FALSE))



setwd("/Users/urszulaczerwinska/Documents/CURIE2015/Data/Biton.Zinovyev/GSE6532/ICAtest/8ics/IC3_BreastBCR_min/ICA_signif_genes")
write.table(All_sign.cl.m.c, file="significative_numerical", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
write.table(row.names(All_sign.cl.m.c), file= "significative_ids", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
write.table(colnames(All_sign.cl.m.c), file= "significative_samples", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)



hist(IC5.IC3.s, col="black")
barplot(IC2.IC3, main="barplot IC1")

length(row.names(IC3.IC3_plus))


IC.IC3_enrichment<-enrichment(row.names(IC1.IC3_min), min_module_size = 10, threshold = 0.05, maps = list(all=cosmos ))
IC.IC3_enrichment[IC.IC3_enrichment$p.value.corrected<0.5,][2:7,]

represent_enrichment(IC.IC3_enrichment[1:6,], plot="bar", scale="log")

IC2.IC3_enrichment<-enrichment(row.names(IC2.IC3_plus), min_module_size = 10, threshold = 0.05, maps = list(all=cosmos))
IC2.IC3_enrichment[,1]
represent_enrichment(IC2.IC3_enrichment, plot="bar", scale="log")



IC3.IC3_enrichment<-enrichment(row.names(IC3.IC3_plus), min_module_size = 5, threshold = 0.5, maps = list(all=cosmos))
IC3.IC3_enrichment[,]


IC4.IC3_enrichment<-enrichment(row.names(IC4.IC3_plus), min_module_size = 5, threshold = 0.05, maps = list(all=cosmos ))
IC4.IC3_enrichment[,]



IC5.IC3_enrichment<-enrichment(row.names(IC5.IC3_plus), min_module_size = 10, threshold = 0.05, maps = list(all=cosmos))
IC5.IC3_enrichment[,1]
represent_enrichment(IC5.IC3_enrichment[1:8,], plot="bar", scale="log")

IC5.IC3_enrichment_m<-enrichment(row.names(IC5.IC3_min), min_module_size = 10, threshold = 0.5, maps = list(Abbas=Abbas, Grigo=Grigo, IRIS=IRIS, Ha=Ha, Palmer=Palmer,Bindea=Bindea,LM22=LM22_gmt))
IC5.IC3_enrichment_m[,]


#sample weights
IC3.samples.names=read.table("/Users/urszulaczerwinska/Documents/CURIE2015/Data/Biton.Zinovyev/GSE6532/ICAtest/8ics/IC3_BreastBCR_min/IC3_min_immune_component_samples", sep="\t", header=F, row.names = NULL)


IC1.IC3.sa<- read.table("/Users/urszulaczerwinska/Documents/CURIE2015/Data/Biton.Zinovyev/GSE6532/ICAtest/8ics/IC3_BreastBCR_min/A_IC3_min_immune_component_numerical_5.num", sep="\t", header=F, row.names = NULL)
IC1.IC3.sa=IC1.IC3.sa[,1:5]
summary(IC1.sa)
row.names(IC1.IC3.sa)=IC3.samples.names[,1]
summary(IC3.samples.names)
head(IC1.IC3.sa)

par(mfrow=c(1,5))
hist(scale(IC1.IC3.sa[,1]))
hist(scale(IC1.IC3.sa[,2]))
hist(scale(IC1.IC3.sa[,3]))
hist(scale(IC1.IC3.sa[,4]))
hist(scale(IC1.IC3.sa[,5]))

barplot(IC1.IC3.sa[,5])


LM22<- read.table("/Users/urszulaczerwinska/Documents/CURIE2015/Data/NewmanCIBERSORT/LM22.txt", sep="\t", header=T, row.names = 1)
row.names(LM22)
sig.rows=intersect(row.names(LM22),row.names(IC1.IC3))
length(sig.rows)
IC.LM22=cbind(IC1.IC3[sig.rows,],IC2.IC3[sig.rows,], IC3.IC3[sig.rows,], IC4.IC3[sig.rows,], IC5.IC3[sig.rows,], LM22[sig.rows,])
IC.LM22.s=scale(IC.LM22)
summary(IC.LM22.s)
IC.LM22.cor=rcorr(IC.LM22.s)

corrplot(IC.LM22.cor$r[1:5,6:27], type="full", order="original", tl.col="black", tl.srt=45, col=colmy(300))
dim(IC.LM22.cor$r)


mLM22.IC3=flattenCorrMatrix(IC.LM22.cor$r, IC.LM22.cor$P)

mLM22.IC3.i=mLM22.IC3[c(11:20,22:26,29:33,37:41,46:50, 56:60,67:71, 79:83, 92:96,106:110,121:125,137:141,154:158,172:176,191:195,211:215, 232:236,254:258,277:281,301:305,326:330),]
mLM22.IC3.i=mLM22.IC3.i[mLM22.IC3.i$p<0.05,]
mLM22.IC3.i=mLM22.IC3.i[abs(mLM22.IC3.i$cor)>0.2,]
mLM22.IC3.o <- mLM22.IC3.i[ order(mLM22.IC3.i$row), ]


#________________________
setwd("/Users/urszulaczerwinska/Documents/CURIE2015/Data/Biton.Zinovyev/GSE6532/ICAtest/8ics/IC3_BreastBCR_min/Immune_sig_gmt")

celTypes=c("B cells naive",	"B cells memory","Plasma cells","T cells CD8","T cells CD4 naive","T cells CD4","memory resting	T cells","CD4 memory activated","T cells follicular helper","T cells regulatory (Tregs)",	"T cells gamma delta","NK cells resting	NK cells activated","Monocytes","Macrophages M0","Macrophages M1",	"Macrophages M2","Dendritic cells resting",	"Dendritic cells activated","Mast cells resting","Mast cells activated", "Eosinophils","Neutrophils")
LM22_bin<- read.table("pre_gmt.txt", sep="\t", header=F, row.names = 1)
for (i in 1:ncol(LM22_bin)){
  print(i)
  #i=3
  print(row.names(LM22_bin[which(LM22_bin[,i]>0),]))
  
  
  write.table(t(c(gsub('([[:punct:]])|\\s+','_', paste("LM22",celTypes[i],sep=" ")), "na",row.names(LM22_bin[which(LM22_bin[,i]>0),]))), file= paste(gsub('([[:punct:]])|\\s+','_',paste("LM22_",celTypes[i],sep="")),".gmt",sep=""), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
  
  
}

#_________________________
CIBERSORT_IC3=read.table("/Users/urszulaczerwinska/Documents/CURIE2015/Data/Biton.Zinovyev/GSE6532/ICAtest/8ics/IC3_BreastBCR_min/CIBERSORT\ deconv/CIBERSORT.IC3_LM22.txt", sep="\t", header=TRUE, row.names = 1)
summary(CIBERSORT_IC3)
head(CIBERSORT_IC3)
CIBERSORT_IC3.o=CIBERSORT_IC3[order(CIBERSORT_IC3$'P_value'),]
CIBERSORT_IC3.o.s=CIBERSORT_IC3.o[CIBERSORT_IC3.o$P_value<0.05,]
dim(CIBERSORT_IC3.o.s)

par(mar=c(7,3,1,1))
    
barplot(apply(CIBERSORT_IC3.o.s[,1:22], 2, FUN=mean), las=2, cex.names=0.5)

panBreast=read.table("/Users/urszulaczerwinska/Documents/CURIE2015/Data/Biton.Zinovyev/GSE6532/ICAtest/8ics/IC3_BreastBCR_min/CIBERSORT\ deconv/panBreast_CIBER.txt", sep="\t", header=F)
colnames(panBreast)=colnames(CIBERSORT_IC3.o.s[,1:22])
barplot(apply(panBreast, 2, FUN=mean), las=2, cex.names=0.5)


CIBERSORT_GSE6532=read.table("/Users/urszulaczerwinska/Documents/CURIE2015/Data/Biton.Zinovyev/GSE6532/CIBERSORT/CIBERSORT.Output_Job3.txt", sep="\t", header=TRUE, row.names = 1)
head(CIBERSORT_GSE6532)
CIBERSORT_GSE6532.s=CIBERSORT_GSE6532[CIBERSORT_GSE6532$P.value<0.05,]
dim(CIBERSORT_GSE6532.s)
head(CIBERSORT_GSE6532.s)
A=apply(CIBERSORT_GSE6532.s[,1:22], 2, FUN=mean)
barplot(as.matrix(cbind(A,A,A)), las=2, cex.names=0.5,col=rainbow(22), legend=TRUE)

#-NMF on the GSE6532
dim(data.tam.t.n.u)
data.tam.t.n.u=data.tam.t.n.u[,2:278]
exp(data.tam.t.n.u)

estim.r <- nmf(exp(data.tam.t.n.u), 2:10, nrun = 10, seed = 123456)
plot(estim.r)
consensusmap(estim.r)

# shuffle original data
V.random <- randomize(exp(data.tam.t.n.u))
# estimate quality measures from the shuffled data (use
# default NMF algorithm)
estim.r.random <- nmf(V.random, 2:10, nrun = 10, seed = 123456)
# plot measures on same graph
plot(estim.r, estim.r.random)


res <- nmf(exp(data.tam.t.n.u), 8, .options = "t")
plot(res)


res.multi.method <- nmf(exp(data.tam.t.n.u), 8, list("brunet", "lee", "ns"), seed = 123456, .options = "t")
plot(res.multi.method)


res.2 <- nmf(exp(data.tam.t.n.u), 8, method="lee", .options='v',seed=123456)

res.3 <- nmf(exp(data.tam.t.n.u), 6, method="lee", .options='v',seed=123456)

H=coef(res.2)
W=basis(res.2)
setwd("/Users/urszulaczerwinska/Documents/CURIE2015/Data/Biton.Zinovyev/GSE6532/NMF")
write.table(coef(res.2), file="samples_coef.txt", quote=FALSE, sep="\t")
write.table(basis(res.2), file="genes_basis.txt", quote=FALSE, sep="\t",col.names = FALSE)
basismap(res.2)
basismap(res.3)
s <- featureScore(res.2)
summary(s)

s <- extractFeatures(res.2)
str(s)
s


NMF.1.s=basis(res.2)[s[[1]],1]
NMF.2.s=basis(res.2)[s[[2]],2]
NMF.3.s=basis(res.2)[s[[3]],3]
NMF.4.s=basis(res.2)[s[[4]],4]
NMF.5.s=basis(res.2)[s[[5]],5]
NMF.6.s=basis(res.2)[s[[6]],6]
NMF.7.s=basis(res.2)[s[[7]],7]
NMF.8.s=basis(res.2)[s[[8]],8]

NMF.1.o=NMF.1.s[order(log(NMF.1.s))]
NMF.2.o=NMF.2.s[order(log(NMF.2.s))]
NMF.3.o=NMF.3.s[order(log(NMF.3.s))]
NMF.4.o=NMF.4.s[order(log(NMF.4.s))]
NMF.5.o=NMF.5.s[order(log(NMF.5.s))]
NMF.6.o=NMF.6.s[order(log(NMF.6.s))]
NMF.7.o=NMF.7.s[order(log(NMF.7.s))]
NMF.8.o=NMF.8.s[order(log(NMF.8.s))]





setwd("/Users/urszulaczerwinska/Documents/CURIE2015/Data/Biton.Zinovyev/GSE6532/NMF/clust")

write.table(NMF.1.o, file="NMF.1_clust.rnk", quote=FALSE, sep="\t",col.names=FALSE)
write.table(NMF.2.o, file="NMF.2_clust.rnk", quote=FALSE, sep="\t",col.names=FALSE)
write.table(NMF.3.o, file="NMF.3_clust.rnk", quote=FALSE, sep="\t",col.names=FALSE)
write.table(NMF.4.o, file="NMF.4_clust.rnk", quote=FALSE, sep="\t",col.names=FALSE)
write.table(NMF.5.o, file="NMF.5_clust.rnk", quote=FALSE, sep="\t",col.names=FALSE)
write.table(NMF.6.o, file="NMF.6_clust.rnk", quote=FALSE, sep="\t",col.names=FALSE)
write.table(NMF.7.o, file="NMF.7_clust.rnk", quote=FALSE, sep="\t",col.names=FALSE)
write.table(NMF.8.o, file="NMF.8_clust.rnk", quote=FALSE, sep="\t",col.names=FALSE)



histhist(log(W[,1]))

barplot(log(W[order(W[,1]),1]))
length(W[order(W[,1]),1])
W[order(W[,1]),1][14699:14709]

names(NMF.1.s)[1:10]


W.1.s=scale(W[,1])
W.1.l.s=scale(log(W[,1]))
hist(W.1.l.s)
W.1.l.t=W.1.l.s[W.1.l.s>0.5,]

W.1.s=scale(W[,8])
W.1.l.s=scale(log(W[,8]))
hist(W.1.l.s)
length(W.1.l.s[W.1.l.s>0.5,])
W.1.l.t=W.1.l.s[W.1.l.s>0.5,]
setwd("/Users/urszulaczerwinska/Documents/CURIE2015/Data/Biton.Zinovyev/GSE6532/NMF/tr.05")

for (i in 1:8){
  
  W.1.l.s=scale(log(W[,i]))
  W.1.l.t=W.1.l.s[W.1.l.s>0.5,]
  W.1.l.o=W.1.l.t[order(W.1.l.t)]
  write.table( W.1.l.o, file=paste("NMF.",i,".rnk",sep=""), quote=FALSE, sep="\t",col.names=FALSE)
  
  
}



W.8.s=scale(W[,8])
W.8.l.s=scale(log(W[,8]))
hist(W.8.l.s)
length(W.8.l.s[W.1.l.s>0.5,])
W.8.l.t=W.1.l.s[W.8.l.s>0.5,]

names(W.8.l.t)
dim(W)
head(data.tam.t.n.u)
NMF.8.data=data.tam.t.n.u[names(W.8.l.t),]

W.6=basis(res.3)
W.l.s=matrix(0, ncol=6, nrow=14709)
for (i in 1:6) {
  i
  W.s=scale(W.6[,i])
  W.l.s[,i]=scale(log(W[,i]))
  


}
summary(W.l.s)


par(mfrow=c(2,3))
setwd("/Users/urszulaczerwinska/Documents/CURIE2015/Data/Biton.Zinovyev/GSE6532/NMF/06.01.16")

#i=3
for (i in 1:6){
#hist(W.l.s[,i])
#length(W.l.s[,i][W.l.s[,i]>0.5])
#W.l.t=W.1.l.s[W.l.s>0.5,]
  W.s=scale(log(W.6[,i]))
  W.l.t=W.s[W.s>0.5,]
  W.l.o=W.l.t[order(W.l.t)]
  print(length(W.l.o))
  #write.table( W.l.o, file=paste("NMF_6_.",i,".rnk",sep=""), quote=FALSE, sep="\t",col.names=FALSE)
  
  
}


#res.multi.method <- nmf(exp(data.tam.t.n.u), 8, list("brunet", "lee", "ns"), seed = 123456, .options = "t")
#plot(res.multi.method)

estim.r <- nmf(exp(NMF.8.data), 2:10, nrun = 10, seed = 123456)
plot(estim.r)

estim.r.2 <- nmf(exp(NMF.8.data), 8:15, nrun = 10, seed = 123456)
plot(estim.r.2)


estim.r.3 <- nmf(exp(NMF.8.data), 14:20, nrun = 10, seed = 123456)
plot(estim.r.3)

estim.r.4 <- nmf(exp(NMF.8.data), 19:24, nrun = 10, seed = 123456)
plot(estim.r.4)

i=4
W.s=scale(log(W.6[,i]))
W.l.t=W.s[W.s>0.5,]
W.l.o=W.l.t[order(W.l.t)]
length(names(W.l.o))
dim(W.6)
head(data.tam.t.n.u)
NMF_6.3.data=data.tam.t.n.u[names(W.l.o),]
NMF_6.4.data=data.tam.t.n.u[names(W.l.o),]

dim(NMF_6.3.data)

estim.r_6.3 <- nmf(exp(NMF_6.3.data), 2:5, nrun = 10, seed = 123456)
plot(estim.r_6.3)

estim.r_6.3_2 <- nmf(exp(NMF_6.3.data), 6:10, nrun = 10, seed = 123456)
plot(estim.r_6.3_2)


res.NMF_6.3 <- nmf(exp(NMF_6.3.data), 8, method="lee", .options='v',seed=123456)
W.6.3_8=basis(res.NMF_6.3)

par(mfrow=c(2,4))
setwd("/Users/urszulaczerwinska/Documents/CURIE2015/Data/Biton.Zinovyev/GSE6532/NMF/07.01.16")

for (i in 1:8){
 W.6.3_8.i=scale(log(W.6.3_8[,i]))
  #hist(W.6.3_8.i)
  print(length(W.6.3_8.i[W.6.3_8.i>0.5]))
  #W.l.t=W.1.l.s[W.l.s>0.5,]
  #W.s=scale(log(W.6[,i]))
  #W.l.t=W.s[W.s>0.5,]
  W.o=W.6.3_8.i[W.6.3_8.i>0.5,, drop=FALSE]
  W.l.o=W.o[order(W.o),,drop=FALSE]
assign(paste("W.6.3_8",i,sep="."),W.l.o)
  #print(length(W.l.o))
  write.table( W.l.o, file=paste("NMF_6.3_8.",i,".rnk",sep=""), quote=FALSE, sep="\t",col.names=FALSE)
  
  
}
for (i in 4:6){
NMF_enrichment<-enrichment(row.names(noquote(paste("W.6.3_8", i, sep="."))), min_module_size = 10, threshold = 0.05, maps = list(all=cosmos))
assign(paste("NMF_enrichment", i, sep="."),NMF_enrichment)
}


NMF_enrichment.1<-enrichment(row.names(W.6.3_8.1), min_module_size = 10, threshold = 0.05, maps = list(all=cosmos))
NMF_enrichment.2<-enrichment(row.names(W.6.3_8.2), min_module_size = 10, threshold = 0.05, maps = list(all=cosmos))
NMF_enrichment.3<-enrichment(row.names(W.6.3_8.3), min_module_size = 10, threshold = 0.05, maps = list(all=cosmos))
NMF_enrichment.4<-enrichment(row.names(W.6.3_8.4), min_module_size = 10, threshold = 0.05, maps = list(all=cosmos))
NMF_enrichment.5<-enrichment(row.names(W.6.3_8.5), min_module_size = 10, threshold = 0.05, maps = list(all=cosmos))
NMF_enrichment.6<-enrichment(row.names(W.6.3_8.6), min_module_size = 10, threshold = 0.05, maps = list(all=cosmos))
NMF_enrichment.7<-enrichment(row.names(W.6.3_8.7), min_module_size = 10, threshold = 0.05, maps = list(all=cosmos))
NMF_enrichment.8<-enrichment(row.names(W.6.3_8.8), min_module_size = 10, threshold = 0.05, maps = list(all=cosmos))


estim.r_6.4 <- nmf(exp(NMF_6.4.data), 2:10, nrun = 10, seed = 123456)
plot(estim.r_6.4)


res.NMF_6.4 <- nmf(exp(NMF_6.4.data), 8, method="lee", .options='v',seed=123456)


res.NMF_6.4 <- nmf(exp(NMF_6.4.data), 8, method="lee", .options='v',seed=123456)
W.6.4_8=basis(res.NMF_6.4)

par(mfrow=c(2,4))
setwd("/Users/urszulaczerwinska/Documents/CURIE2015/Data/Biton.Zinovyev/GSE6532/NMF/07.01.16")

for (i in 1:8){
  W.6.4_8.i=scale(log(W.6.4_8[,i]))
 # hist(W.6.4_8.i)
  #print(length(W.6.4_8.i[W.6.4_8.i>0.5]))
 # W.l.t=W.1.l.s[W.l.s>0.5,]
  #W.s=scale(log(W.6[,i]))
  #W.l.t=W.s[W.s>0.5,]
  W.o=W.6.4_8.i[W.6.4_8.i>0.5,, drop=FALSE]
  W.l.o=W.o[order(W.o),,drop=FALSE]
  assign(paste("W.6.4_8",i,sep="."),W.l.o)
  #print(length(W.l.o))
  write.table( W.l.o, file=paste("NMF_6.4_8.",i,".rnk",sep=""), quote=FALSE, sep="\t",col.names=FALSE)
  
  
}

NMF_enrichment_4.1<-enrichment(row.names(W.6.4_8.1), min_module_size = 10, threshold = 0.05, maps = list(all=cosmos))
NMF_enrichment_4.2<-enrichment(row.names(W.6.4_8.2), min_module_size = 10, threshold = 0.05, maps = list(all=cosmos))
NMF_enrichment_4.3<-enrichment(row.names(W.6.4_8.3), min_module_size = 10, threshold = 0.05, maps = list(all=cosmos))
NMF_enrichment_4.4<-enrichment(row.names(W.6.4_8.4), min_module_size = 10, threshold = 0.05, maps = list(all=cosmos))
NMF_enrichment_4.5<-enrichment(row.names(W.6.4_8.5), min_module_size = 10, threshold = 0.05, maps = list(all=cosmos))
NMF_enrichment_4.6<-enrichment(row.names(W.6.4_8.6), min_module_size = 10, threshold = 0.05, maps = list(all=cosmos))
NMF_enrichment_4.7<-enrichment(row.names(W.6.4_8.7), min_module_size = 10, threshold = 0.05, maps = list(all=cosmos))
NMF_enrichment_4.8<-enrichment(row.names(W.6.4_8.8), min_module_size = 10, threshold = 0.05, maps = list(all=cosmos))

##############2 IC ICA #############################

setwd("/Users/urszulaczerwinska/Documents/CURIE2015/Data/Biton.Zinovyev/GSE6532/ICAtest/8ics/IC3_BreastBCR_min")

IC3_2<- read.table("./S_IC3_min_immune_component_numerical_2.annot", sep="\t", header=T, row.names = 1)
summary(IC3_2)
IC1.IC3_2=IC3_2[,1,drop=FALSE]
IC2.IC3_2=IC3_2[,2,drop=FALSE]
row.names(IC2.IC3_2)
names(IC2.IC3_2)="IC2"

IC1.IC3_2.o=IC1.IC3_2[order(IC1.IC3_2),,drop=FALSE]
IC2.IC3_2.o=IC2.IC3_2[order(IC2.IC3_2),,drop=FALSE]
IC1.IC3_2.o.s=scale(IC1.IC3_2.o)
IC2.IC3_2.o.s=scale(IC2.IC3_2.o)
par(mfrow=c(1,2))
hist(IC1.IC3_2.o.s)
hist(IC2.IC3_2.o.s)

IC1.IC3_2.min=IC1.IC3_2.o.s[IC1.IC3_2.o.s< -2,,drop=FALSE]
IC1.IC3_2.plus=IC1.IC3_2.o.s[IC1.IC3_2.o.s>2,,drop=FALSE]

length(IC1.IC3_2.min)
length(IC1.IC3_2.plus)



IC2.IC3_2.min=IC2.IC3_2.o.s[IC2.IC3_2.o.s< -2,,drop=FALSE]
IC2.IC3_2.plus=IC2.IC3_2.o.s[IC2.IC3_2.o.s>2,,drop=FALSE]

length(IC2.IC3_2.min)
length(IC2.IC3_2.plus)



IC1.IC3_2_enrichment_min<-enrichment(row.names(IC1.IC3_2.min), min_module_size = 10, threshold = 0.05, maps = list(all=cosmos))
IC1.IC3_2_enrichment_plus<-enrichment(row.names(IC1.IC3_2.plus), min_module_size = 10, threshold = 0.05, maps = list(all=cosmos))

IC2.IC3_2_enrichment_min<-enrichment(row.names(IC2.IC3_2.min), min_module_size = 10, threshold = 0.05, maps = list(all=cosmos))
IC2.IC3_2_enrichment_plus<-enrichment(row.names(IC2.IC3_2.plus), min_module_size = 10, threshold = 0.05, maps = list(all=cosmos))

sig.rows=intersect(row.names(LM22),row.names(IC1.IC3_2.o.s))
length(sig.rows)
IC_2.LM22=cbind(IC1.IC3_2.o.s[sig.rows,],IC2.IC3_2.o.s[sig.rows,],LM22[sig.rows,])
IC_2.LM22.s=scale(IC_2.LM22)
summary(IC_2.LM22.s)
IC_2.LM22.cor=rcorr(IC_2.LM22.s)

corrplot(IC_2.LM22.cor$r[1:2,3:24], type="full", order="original", tl.col="black", tl.srt=45, col=colmy(300))
dim(IC.LM22.cor$r)


setwd("/Users/urszulaczerwinska/Documents/CURIE2015/Data/Biton.Zinovyev/GSE6532/ICAtest/8ics/IC3_BreastBCR_min")

IC3_6<- read.table("./S_IC3_min_immune_component_numerical_6.annot", sep="\t", header=T, row.names = 1)
summary(IC3_6)
IC3_6=IC3_6[,1:6]
names(IC3_6)=c("IC1","IC2", "IC3", "IC4", "IC5", "IC6")

IC3_6.s=scale(IC3_6)

par(mfrow=c(2,3))

for (i in 1:6)
{
  hist(IC3_6.s[,i])
  
}


IC_6.LM22=cbind(IC3_6.s[sig.rows,],LM22[sig.rows,])
IC_6.LM22.s=scale(IC_6.LM22)
summary(IC_6.LM22.s)
IC_6.LM22.cor=rcorr(IC_6.LM22.s)

corrplot(IC_6.LM22.cor$r[1:6,7:28], type="full", order="original", tl.col="black", tl.srt=45, col=colmy(300))
dim(IC.LM22.cor$r)


IC1.IC3_6.min=IC3_6.s[IC3_6.s[,1]< -2,1,drop=FALSE]
IC1.IC3_6_enrichment_min<-enrichment(row.names(IC1.IC3_6.min), min_module_size = 10, threshold = 0.05, maps = list(all=cosmos))

IC2.IC3_6.min=IC3_6.s[IC3_6.s[,2]< -2,2,drop=FALSE]
IC2.IC3_6_enrichment_min<-enrichment(row.names(IC2.IC3_6.min), min_module_size = 10, threshold = 0.05, maps = list(all=cosmos))

IC3.IC3_6.min=IC3_6.s[IC3_6.s[,3]< -2,3,drop=FALSE]
IC3.IC3_6_enrichment_min<-enrichment(row.names(IC3.IC3_6.min), min_module_size = 10, threshold = 0.05, maps = list(all=cosmos))


IC3.IC3_6.plus=IC3_6.s[IC3_6.s[,3] > 2, 3, drop=FALSE]
IC3.IC3_6_enrichment_plus<-enrichment(row.names(IC3.IC3_6.plus), min_module_size = 10, threshold = 0.05, maps = list(all=cosmos))


IC4.IC3_6.min=IC3_6.s[IC3_6.s[,4]< -2,3,drop=FALSE]
IC4.IC3_6_enrichment_min<-enrichment(row.names(IC4.IC3_6.min), min_module_size = 10, threshold = 0.05, maps = list(all=cosmos))

IC5.IC3_6.plus=IC3_6.s[IC3_6.s[,5] > 2, 5, drop=FALSE]
IC5.IC3_6_enrichment_plus<-enrichment(row.names(IC5.IC3_6.plus), min_module_size = 10, threshold = 0.05, maps = list(all=cosmos))

IC6.IC3_6.plus=IC3_6.s[IC3_6.s[,6] > 2, 6, drop=FALSE]
IC6.IC3_6_enrichment_plus<-enrichment(row.names(IC6.IC3_6.plus), min_module_size = 10, threshold = 0.05, maps = list(all=cosmos))




An_IC1=annot.tam[match(row.names(IC1.IC3_6.min), annot.tam[,7], nomatch = 0), c(7,10)]
An_IC2=annot.tam[match(row.names(IC2.IC3_6.min), annot.tam[,7], nomatch = 0), c(7,10)]
An_IC3=annot.tam[match(row.names(IC3.IC3_6.min), annot.tam[,7], nomatch = 0), c(7,10)]
An_IC4=annot.tam[match(row.names(IC4.IC3_6.min), annot.tam[,7], nomatch = 0), c(7,10)]
An_IC5=annot.tam[match(row.names(IC5.IC3_6.plus), annot.tam[,7], nomatch = 0), c(7,10)]
An_IC6=annot.tam[match(row.names(IC6.IC3_6.plus), annot.tam[,7], nomatch = 0), c(7,10)]





################# MineICA ##########################
library(Biobase)
library(plyr)
library(ggplot2)
library(foreach)
library(xtable)
library(biomaRt)
library(GOstats)
library(cluster)
library(marray)
library(mclust)
library(RColorBrewer)
library(igraph)
library(Rgraphviz)
library(graph)
library(colorspace)
library(annotate)
library(scales)
library(gtools)

library(MineICA)

## load Mainz expression data and sample annotations.
library(breastCancerMAINZ)
data(mainz)
show(mainz)
mainz <- selectFeatures_IQR(mainz,10000)
library(JADE)
## Features are mean-centered before ICA computation
exprs(mainz) <- t(apply(exprs(mainz),1,scale,scale=FALSE))
colnames(exprs(mainz)) <- sampleNames(mainz)
## run ICA-JADE
resJade <- runICA(X=exprs(mainz), nbComp=5, method = "JADE", maxit=10000)

library(fastICA)
res <- clusterFastICARuns(X=exprs(mainz), nbComp=5, alg.type="deflation", nbIt=10, funClus="hclust", method="average")

str(res)
## build params
setwd("/Users/urszulaczerwinska/Documents/CURIE2015/MineICA")
params <- buildMineICAParams(resPath="./", selCutoff=3, pvalCutoff=0.05)

## load annotation package
library(hgu133a.db)

ls("package:hgu133a.db")
mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
ensembl = useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host = "www.ensembl.org")


listFilters(ensembl)[,]
listAttributes(ensembl)[grep(x=listAttributes(mart)[,1],pattern="affy")[1:5],]

## Define typeID, Mainz data originate from affymetrix HG-U133a microarray
## and are indexed by probe sets.
## The probe sets are annotated into Gene Symbols
typeIDmainz <- c(geneID_annotation="SYMBOL", geneID_biomart="hgnc_symbol", featureID_biomart=c("affy_hg_u133a","affy_hg_u133b","affy_hg_u133_plus_2"))
## define the reference samples if any, here no normal sample is available
refSamplesMainz <- character(0)
resBuild <- buildIcaSet(params=params, A=data.frame(resJade$A), S=data.frame(resJade$S), dat=exprs(mainz), pData=pData(mainz), refSamples=refSamplesMainz,annotation="hgu133a.db", typeID= typeIDmainz, chipManu = "affymetrix", mart=ensembl)
icaSetMainz <- resBuild$icaSet
str(icaSetMainz)
params <- resBuild$params
annot <- pData(icaSetMainz)
varLabels(icaSetMainz)[1:5]
icaSetMainz$grade[1:5]
featureNames(icaSetMainz)[1:5] # probe set ids
geneNames(icaSetMainz)[1:5] #gene symbols
sampleNames(icaSetMainz)[1:5]
head(dat(icaSetMainz))
head(datByGene(icaSetMainz))
A(icaSetMainz)
S(icaSetMainz)
SByGene(icaSetMainz)
nbComp(icaSetMainz)
compNames(icaSetMainz)
indComp(icaSetMainz)

hist(S(icaSetMainz)[,1])
## Extract the contributing genes
contrib <- selectContrib(icaSetMainz, cutoff=3, level="genes")
## Show the first contributing genes of the first and third components
sort(abs(contrib[[1]]),decreasing=TRUE)[1:10]
sort(abs(contrib[[3]]),decreasing=TRUE)[1:10]

## One can also want to apply different cutoffs depending on the components
## for example using the first 4 components:
contrib <- selectContrib(icaSetMainz[,,1:4], cutoff=c(4,4,4,3), level="genes")

## extract sample contributions and gene projections of the second component
comp2 <- getComp(icaSetMainz, level="genes", ind=2)
## access the sample contributions
comp2$contrib[1:5]

## access the gene projections
comp2$proj[1:5]

## select the annotations of interest
varLabels(icaSetMainz)

# restrict the phenotype data to the variables of interest
keepVar <- c("age","er","grade")
# specify the variables that should be treated as character
icaSetMainz$er <- c("0"="ER-","1"="ER+")[as.character(icaSetMainz$er)]
icaSetMainz$grade <- as.character(icaSetMainz$grade)
## Run the analysis of the ICA decomposition
# only enrichment in KEGG gene sets are tested
runAn(params=params, icaSet=icaSetMainz, writeGenesByComp = TRUE,mart = ensembl, keepVar=keepVar, dbGOstats = "KEGG")


resW <- writeProjByComp(icaSet=icaSetMainz, params=params, mart=ensembl, level='genes', selCutoffWrite=2.5)
str(resW)
head(resW$listAnnotComp[[1]])

head(resW$nbOccInComp)

#Heatmaps
## selection of the variables we want to display on the heatmap
keepVar <- c("er","grade")
## For the second component, select contributing genes using a threshold of 3
## on the absolute projection values,
## heatmap with dendrogram
resH <- plot_heatmapsOnSel(icaSet = icaSetMainz, selCutoff = 3, level = "genes", keepVar = keepVar, doSamplesDendro = TRUE, doGenesDendro = TRUE, keepComp = 2,heatmapCol = maPalette(low = "blue", high = "red", mid = "yellow", k=44), file = "heatmapWithDendro", annot2col=annot2col(params))
## heatmap where genes and samples are ordered by contribution values
resH <- plot_heatmapsOnSel(icaSet = icaSetMainz, selCutoff = 3, level = "genes",keepVar = keepVar,doSamplesDendro = FALSE, doGenesDendro = FALSE, keepComp = 2, heatmapCol = maPalette(low = "blue", high = "red", mid = "yellow", k=44), file = "heatmapWithoutDendro", annot2col=annot2col(params))

### Qualitative variables
## Compute Wilcoxon and Kruskall-Wallis tests to compare the distribution
## of the samples according to their grade and ER status on all components.
resQual <- qualVarAnalysis(params=params, icaSet=icaSetMainz, keepVar=c("er","grade"), adjustBy="none", typePlot="boxplot", path="qualVarAnalysis/", filename="qualVar")

### Quantitative variables
## Compute pearson correlations between variable ✬age✬ and the sample contributions
## on all components.
## We are interested in correlations exceeding 0.3 in absolute value, and plots will only be drawn
## for correlations exceeding this threshold.
resQuant <- quantVarAnalysis(params=params, icaSet=icaSetMainz, keepVar="age", typeCor="pearson", cutoffOn="cor", cutoff=0.3, adjustBy="none", path="quantVarAnalysis/", filename="quantVar", typeImage = "svg")

resmix <- plotAllMix(A=A(icaSetMainz), nbMix=2, nbBreaks=50)

## plot the positions of the samples on the second component according to their ER status
## (in a file "er.pdf")
plotPosAnnotInComp(icaSet=icaSetMainz, params=params, keepVar=c("er"), keepComp=2,funClus="Mclust")


## clustering of the samples in 1-dim using the vector
## of sample contributions of the two first components
## and Gaussian mixture modeling (Mclust)
clus1 <- clusterSamplesByComp(params=params, icaSet=icaSetMainz[,,,1:2],funClus="Mclust", clusterOn="A", nbClus=2, filename="comp1Mclust")

## The obtained clusters are written in the file "comp1Mclus.txt" of the result path.
clus1$clus[[2]][1:5]

clus2 <- clusterSamplesByComp_multiple(params=params, icaSet=icaSetMainz[,,1:2],funClus="kmeans", clusterOn=c("A","S"), level="features",nbClus=2, filename="comparKmeans")
## The obtained clusters and their comparison with adjusted Rand indices are written
## in file "comparKmeans.txt" of the result path.
## Both clustering results are stored in a common data.frame
head(clus2$clus)

## Access Rand index
clus2$comparClus


## Test the association between the clustering obtained by Mclust for the first
## component and the variables:
clus2var <- clusVarAnalysis(icaSet=icaSetMainz[,,1:2], params=params, keepVar=c("er","grade"),resClus=clus1$clus, funClus="Mclust", adjustBy="none",doPlot=TRUE, path="clus2var/", filename="resChitests-Mcluscomp1")
## Look at the filename which contains p-values and links to the barplots
## p-values are also contained in the ouput of the function:
clus2var


########BreastBCR

## build params
setwd("/Users/urszulaczerwinska/Documents/CURIE2015/Data/Biton.Zinovyev/GSE6532/MineICA_BreastBCR")
params <- buildMineICAParams(resPath="./", selCutoff=2, pvalCutoff=0.05)

refSamplesMainz <- character(0)
#SLOWWWW AND ERROR res_breast_bcr_8<- clusterFastICARuns(X=data.tam.t, nbComp=8, alg.type="deflation", nbIt=100, funClus="hclust", method="average")


as.character(data.frame(annot.tam.o)[,6])
annot.tam.uni=na.omit(subset(annot.tam.o[,6], !duplicated(annot.tam.o[,6])))
myrows_ann.2=intersect(names(annot.tam.uni),row.names(S.my))

             
A.raw=read.table("./A_Breast_BCR.numerical.txt_8.num", header=FALSE, sep="\t", row.names=NULL)




A.raw.2=A.raw[,-9]
dim(A.raw.2)
row.names(A.raw.2)=t(colnames(data.tam.t.o))

S.raw=read.table("./S_Breast_BCR.numerical.txt_8.num", header=FALSE, sep="\t", row.names=NULL)
S.raw.2=S.raw[,-9]
dim(S.raw.2)
row.names(S.raw.2)=row.names(data.tam.t.o)
myrows_ann.2=intersect(names(annot.tam.uni),row.names(S.my))

S.raw.3=data.frame(S.raw.2[myrows_ann.2,])
row.names(S.raw.3)=annot.tam.uni
data.tam.t.o.d=data.frame(data.tam.t.o)
dat.raw=data.frame(data.tam.t.o.d[myrows_ann.2,])
row.names(dat.raw)=annot.tam.uni



A.my=data.frame(A.raw.2)
S.my=data.frame(S.raw.2[myrows_ann.2,])
resBreast <- buildIcaSet(params=params, A=A.my, S=S.raw.3, dat=dat.raw, pData=data.frame(demo.tam), refSamples=refSamplesMainz, alreadyAnnot=TRUE)
icaSetBreast <- resBreast$icaSet
params <- resBreast$params

featureNames(icaSetBreast)[1:5] # probe set ids
geneNames(icaSetBreast)
head(dat(icaSetBreast))
head(datByGene(icaSetBreast))
A(icaSetBreast)
S(icaSetBreast)
SByGene(icaSetBreast)
nbComp(icaSetBreast)
compNames(icaSetBreast)
indComp(icaSetBreast)

hist(S(icaSetBreast)[,1])
## Extract the contributing genes
contrib <- selectContrib(icaSetBreast, cutoff=3, level="genes")
## Show the first contributing genes of the first and third components
sort(abs(contrib[[1]]),decreasing=TRUE)[1:10]
sort(abs(contrib[[3]]),decreasing=TRUE)[1:10]

icaSetBreast$er <- c("0"="ER-","1"="ER+")[as.character(icaSetBreast$er)]
icaSetBreast$grade <- as.character(icaSetBreast$grade)
icaSetBreast$age
icaSetBreast$node <- as.character(icaSetBreast$node)
icaSetBreast$series <- as.character(icaSetBreast$series)


runAn(params=params, icaSet=icaSetBreast, writeGenesByComp = TRUE,mart = ensembl, keepVar=c("series","age","grade","size","er","pgr","node","t.rfs","e.rfs","t.dmf","e.dmfs"), dbGOstats = "KEGG")
resQuant <- quantVarAnalysis(params=params, icaSet=icaSetBreast, keepVar="age", typeCor="pearson", cutoffOn="cor", cutoff=0.2, adjustBy="none", path="quantVarAnalysis/", filename="quantVar", typeImage = "svg")



############Immune_IC3min
setwd("/Users/urszulaczerwinska/Documents/CURIE2015/Data/Biton.Zinovyev/GSE6532/ICAtest/8ics/IC3_BreastBCR_min/")
W_IC3=read.table(file="./IC3_min_immune_component_numerical", sep="\t", header=FALSE, row.names=NULL)
summary(W_IC3)
dim(W_IC3)
str(W_IC3)
colnames(W_IC3)=row.names(A.raw.2)

kegg_IC3=read.table(file="./IC3_min_immune_component_hugo_ids.txt", sep="\t",  header=FALSE, row.names=NULL)


row.names(W_IC3)=kegg_IC3[,1]
summary
W_IC3=data.frame(W_IC3)

A.raw_IC3=read.table(file="./A_IC3_min_immune_component_numerical_5.num", sep="\t",  header=FALSE, row.names=NULL)
summary(A.raw_IC3)
A.raw_IC3=A.raw_IC3[,-6]
row.names(A.raw_IC3)=row.names(A.raw.2)
A_IC3=data.frame(A.raw_IC3)


S.raw_IC3=read.table(file="./S_IC3_min_immune_component_numerical_5.num", sep="\t",  header=FALSE, row.names=NULL)
summary(S.raw_IC3)
S.raw_IC3=S.raw_IC3[,-6]
row.names(S.raw_IC3)=kegg_IC3[,1]
S_IC3=data.frame(S.raw_IC3)


## build params
setwd("/Users/urszulaczerwinska/Documents/CURIE2015/Data/Biton.Zinovyev/GSE6532/MineICA_BreastBCR_IC3_immune")
params <- buildMineICAParams(resPath="./", selCutoff=2, pvalCutoff=0.05)

refSamplesMainz <- character(0)
keepVar==c("series","age","grade","size","er","pgr","node","t.rfs","e.rfs","t.dmf","e.dmfs")
resBreast_IC3 <- buildIcaSet(params=params, A=A_IC3, S=S_IC3, dat=W_IC3, pData=data.frame(demo.tam), refSamples=refSamplesMainz, alreadyAnnot=TRUE)
icaSetBreast_IC3 <- resBreast_IC3$icaSet
params <- resBreast_IC3$params

featureNames(icaSetBreast_IC3)[1:5] # probe set ids
geneNames(icaSetBreast_IC3)


icaSetBreast_IC3$er <- c("0"="ER-","1"="ER+")[as.character(icaSetBreast_IC3$er)]
icaSetBreast_IC3$grade <- as.character(icaSetBreast_IC3$grade)
icaSetBreast_IC3$age
icaSetBreast_IC3$node <- c("0"="1","1"="2")[as.character(icaSetBreast_IC3$node)]
icaSetBreast_IC3$series <- as.character(icaSetBreast_IC3$series)

runAn(params=params, icaSet=icaSetBreast_IC3, writeGenesByComp = TRUE,mart = ensembl, keepVar=c("series","age","grade","size","er","pgr","node","t.rfs","e.rfs","t.dmf","e.dmfs"), dbGOstats = "KEGG")
resQuant <- quantVarAnalysis(params=params, icaSet=icaSetBreast_IC3, keepVar="age", typeCor="pearson", cutoffOn="cor", cutoff=0.2, adjustBy="none", path="quantVarAnalysis/", filename="quantVar", typeImage = "svg")

plot(scale(S_IC3[,1]),scale(S_IC3[,5]),col=icaSetBreast_IC3$er)

library(ggplot2)
plot1=qplot(scale(A_IC3[,4]), scale(A_IC3[,5]), data=A_IC3, colour=as.character(icaSetBreast_IC3$grade), size=I(10), alpha=I(0.5))
plot1+theme_bw() +
  
  #eliminates background, gridlines, and chart border
  theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_blank()
  ) +
 theme(axis.line = element_line(color = 'black'))

  
p <- ggplot(A_IC3, aes(icaSetBreast_IC3$grade, scale(A_IC3[,5])))+ geom_boxplot() + geom_jitter(colour=icaSetBreast_IC3$grade)
p

p3 <- ggplot(A_IC3, aes(icaSetBreast_IC3$t.dmf, scale(A_IC3[,2])))+ geom_point() 
p3

cor.test(icaSetBreast_IC3$e.dmfs, scale(A_IC3[,2]), type="spearman")

demo.tam$grade=as.factor(demo.tam$grade)
demo.tam$er=as.factor(demo.tam$er)
demo.tam$e.dmfs=as.factor(demo.tam$e.dmfs)


summary(demo.tam)


##kknn
library(ggplot2)

df     <- data.frame(iris)                   # iris dataset
pca    <- prcomp(df[,1:4], retx=T, scale.=T) # scaled pca [exclude species col]
scores <- pca$x[,1:3]                        # scores for first three PC's

# k-means clustering [assume 3 clusters]
km     <- kmeans(A.sel, centers=3, nstart=5)
ggdata <- data.frame(A.sel, Cluster=km$cluster, Grade=A.sel.an$grade)

# stat_ellipse is not part of the base ggplot package
source("https://raw.github.com/low-decarie/FAAV/master/r/stat-ellipse.R") 

ggplot(ggdata) +
  geom_point(aes(x=scale(A.sel.an[,1]), y=scale(A.sel.an[,2]), color=factor(Grade)), size=I(10), alpha=I(0.6)) +
  stat_ellipse(aes(x=scale(A.sel.an[,1]),y=scale(A.sel.an[,2]),fill=factor(Grade)),
               geom="polygon", level=0.95, alpha=0.15) +
  guides(color=guide_legend("Grade"),fill=guide_legend("Grade")) + 
  theme_bw() +
  
  #eliminates background, gridlines, and chart border
  theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_blank()
  ) +
  theme(axis.line = element_line(color = 'black'))


A.sel=A(icaSetBreast_IC3)[,c(3,5)]
A.sel.an=cbind(A.sel,icaSetBreast_IC3$grade, icaSetBreast_IC3$age)
summary(A.sel.an)
names(A.sel.an)[3]="grade"
names(A.sel.an)[4]="age"


m <- dim(A.sel.an)[1]
val <- sample(1:m, size = round(m/3), replace = FALSE,  prob = rep(1/m, m)) 
A.se.an.learn <- A.sel.an[-val,]
A.se.an.valid <- A.sel.an[val,]

A.4.5.kknn <- kknn(grade~., A.se.an.learn, A.se.an.valid, distance = 1, kernel = "optimal")
summary(A.4.5.kknn)
fit <- fitted(A.4.5.kknn)
table(A.se.an.valid$grade, fit)
pcol <- as.character(as.numeric(A.se.an.valid$grade))
pairs(A.se.an.valid[1:2], pch = pcol, col = c("green3", "red")
      [(A.se.an.valid$grade != fit)+1])




#######ROMA

setwd("/Users/urszulaczerwinska/Documents/CURIE2015/Data/Biton.Zinovyev/GSE6532/ICAtest/8ics/IC3_BreastBCR_min/ROMA_test_01.02.16/resROMA.02.02.16_ICs")
projIC_imm=read.table(file="./moduletable_simple.txt", sep="\t", header=TRUE, row.names=1)
projIC_imm=projIC_imm[,1:277]

summary(projIC_imm)
dim(projIC_imm)
str(projIC_imm)

ggdata <- data.frame(t(projIC_imm), Grade=A.sel.an$grade, Age=A.sel.an$age)

ggplot(ggdata) +
  geom_point(aes(x=scale(IC5), y=scale(IC1), color=factor(Grade)), size=I(10), alpha=I(0.6)) +
  stat_ellipse(aes(x=scale(IC5),y=scale(IC1),fill=factor(Grade)),
               geom="polygon", level=0.95, alpha=0.15) +
  guides(color=guide_legend("Grade"),fill=guide_legend("Grade")) + 
  theme_bw() +
  
  #eliminates background, gridlines, and chart border
  theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_blank()
  ) +
  theme(axis.line = element_line(color = 'black'))

boxplot(t(projIC_imm[1,]),t(projIC_imm[2,]))

X11()
boxplot(grade ~ IC5, IC3, data=ggdata)

A_IC3.nona=na.omit(data.frame(A_IC3,icaSetBreast_IC3$grade))
p <- ggplot(A_IC3.nona, aes(icaSetBreast_IC3.grade, scale(A_IC3.nona[,5])))+ geom_boxplot(aes(fill=icaSetBreast_IC3.grade),alpha=I(0.6)) + geom_jitter(colour="black",size=I(4), alpha=I(0.6))
p=p + 
  
  labs(title=paste("Grade by Tcell - ICA onlycorr\n p.val=",round(k1$p.value,3)), x="Grade", y="T-cell") +
  theme_bw() +
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  #eliminates background, gridlines, and chart border
  theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_blank()
  ) +
  theme(axis.line = element_line(color = 'black'))


ggdata=na.omit(ggdata)
p5 <- ggplot( ggdata, aes(Grade, -IC5))+ geom_boxplot(aes(fill=Grade),alpha=I(0.6)) + geom_jitter(colour="black",size=I(4), alpha=I(0.6))
p5=p5 +
  labs(title=paste("Grade by Tcell - ROMA corr\n p.val=",round(k2$p.value,3)), x="Grade", y="T-cell") +
  theme_bw() +
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  #eliminates background, gridlines, and chart border
  theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_blank()
  ) +
  theme(axis.line = element_line(color = 'black'))

grid.arrange(p, p5, ncol=2)

k1=kruskal.test(A_IC3.nona[,5] ~ icaSetBreast_IC3.grade, data = A_IC3.nona) 
k2=kruskal.test(IC5 ~ Grade, data = ggdata) 

round(k1$p.value,3)

library(FSA)
dunnTest(A_IC3.nona[,5] ~ icaSetBreast_IC3.grade, data = A_IC3.nona, method="fdr")




k1=kruskal.test(A_IC3.nona[,1] ~ icaSetBreast_IC3.grade, data = A_IC3.nona) 
k2=kruskal.test(IC1 ~ Grade, data = ggdata) 
A_IC3.nona=na.omit(data.frame(A_IC3,icaSetBreast_IC3$grade))
p2 <- ggplot(A_IC3.nona, aes(icaSetBreast_IC3.grade, scale(A_IC3.nona[,1])))+ geom_boxplot(aes(fill=icaSetBreast_IC3.grade),alpha=I(0.6)) + geom_jitter(colour="black",size=I(4), alpha=I(0.6))
p2=p2 + 
  
  labs(title=paste("Grade by IC1 - ICA onlycorr\n p.val=",round(k1$p.value,3)), x="Grade", y="T-cell") +
  theme_bw() +
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  #eliminates background, gridlines, and chart border
  theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_blank()
  ) +
  theme(axis.line = element_line(color = 'black'))


ggdata=na.omit(ggdata)
p2b <- ggplot( ggdata, aes(Grade, -IC1))+ geom_boxplot(aes(fill=Grade),alpha=I(0.6)) + geom_jitter(colour="black",size=I(4), alpha=I(0.6))
p2b=p2b +
  labs(title=paste("Grade by IC1 - ROMA corr\n p.val=",round(k2$p.value,3)), x="Grade", y="T-cell") +
  theme_bw() +
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  #eliminates background, gridlines, and chart border
  theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_blank()
  ) +
  theme(axis.line = element_line(color = 'black'))

grid.arrange(p2, p2b, ncol=2)

k1=kruskal.test(A_IC3.nona[,1] ~ icaSetBreast_IC3.grade, data = A_IC3.nona) 
k2=kruskal.test(IC1 ~ Grade, data = ggdata) 
A_IC3.nona=na.omit(data.frame(A_IC3,icaSetBreast_IC3$grade))
p2 <- ggplot(A_IC3.nona, aes(icaSetBreast_IC3.grade, scale(A_IC3.nona[,1])))+ geom_boxplot(aes(fill=icaSetBreast_IC3.grade),alpha=I(0.6)) + geom_jitter(colour="black",size=I(4), alpha=I(0.6))
p2=p2 + 
  
  labs(title=paste("Grade by IC1 - ICA onlycorr\n p.val=",round(k1$p.value,3)), x="Grade", y="T-cell") +
  theme_bw() +
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  #eliminates background, gridlines, and chart border
  theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_blank()
  ) +
  theme(axis.line = element_line(color = 'black'))


ggdata=na.omit(ggdata)
p2b <- ggplot( ggdata, aes(Grade, -IC1))+ geom_boxplot(aes(fill=Grade),alpha=I(0.6)) + geom_jitter(colour="black",size=I(4), alpha=I(0.6))
p2b=p2b +
  labs(title=paste("Grade by IC1 - ROMA corr\n p.val=",round(k2$p.value,3)), x="Grade", y="T-cell") +
  theme_bw() +
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  #eliminates background, gridlines, and chart border
  theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_blank()
  ) +
  theme(axis.line = element_line(color = 'black'))

grid.arrange(p2, p2b, ncol=2)

k1=kruskal.test(A_IC3.nona[,1] ~ icaSetBreast_IC3.grade, data = A_IC3.nona) 
k2=kruskal.test(IC1 ~ Grade, data = ggdata) 
A_IC3.nona=na.omit(data.frame(A_IC3,icaSetBreast_IC3$grade))
p2 <- ggplot(A_IC3.nona, aes(icaSetBreast_IC3.grade, scale(A_IC3.nona[,1])))+ geom_boxplot(aes(fill=icaSetBreast_IC3.grade),alpha=I(0.6)) + geom_jitter(colour="black",size=I(4), alpha=I(0.6))
p2=p2 + 
  
  labs(title=paste("Grade by IC1 - ICA onlycorr\n p.val=",round(k1$p.value,3)), x="Grade", y="T-cell") +
  theme_bw() +
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  #eliminates background, gridlines, and chart border
  theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_blank()
  ) +
  theme(axis.line = element_line(color = 'black'))


ggdata=na.omit(ggdata)
p2b <- ggplot( ggdata, aes(Grade, -IC1))+ geom_boxplot(aes(fill=Grade),alpha=I(0.6)) + geom_jitter(colour="black",size=I(4), alpha=I(0.6))
p2b=p2b +
  labs(title=paste("Grade by IC1 - ROMA corr\n p.val=",round(k2$p.value,3)), x="Grade", y="T-cell") +
  theme_bw() +
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  #eliminates background, gridlines, and chart border
  theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_blank()
  ) +
  theme(axis.line = element_line(color = 'black'))

grid.arrange(p2, p2b, ncol=2)

k1=kruskal.test(A_IC3.nona[,1] ~ icaSetBreast_IC3.grade, data = A_IC3.nona) 
k2=kruskal.test(IC1 ~ Grade, data = ggdata) 
A_IC3.nona=na.omit(data.frame(A_IC3,icaSetBreast_IC3$grade))
p2 <- ggplot(A_IC3.nona, aes(icaSetBreast_IC3.grade, scale(A_IC3.nona[,1])))+ geom_boxplot(aes(fill=icaSetBreast_IC3.grade),alpha=I(0.6)) + geom_jitter(colour="black",size=I(4), alpha=I(0.6))
p2=p2 + 
  
  labs(title=paste("Grade by IC1 - ICA onlycorr\n p.val=",round(k1$p.value,3)), x="Grade", y="T-cell") +
  theme_bw() +
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  #eliminates background, gridlines, and chart border
  theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_blank()
  ) +
  theme(axis.line = element_line(color = 'black'))


ggdata=na.omit(ggdata)
p2b <- ggplot( ggdata, aes(Grade, -IC1))+ geom_boxplot(aes(fill=Grade),alpha=I(0.6)) + geom_jitter(colour="black",size=I(4), alpha=I(0.6))
p2b=p2b +
  labs(title=paste("Grade by IC1 - ROMA corr\n p.val=",round(k2$p.value,3)), x="Grade", y="T-cell") +
  theme_bw() +
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  #eliminates background, gridlines, and chart border
  theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_blank()
  ) +
  theme(axis.line = element_line(color = 'black'))

grid.arrange(p2, p2b, ncol=2)

k1=kruskal.test(A_IC3.nona[,1] ~ icaSetBreast_IC3.grade, data = A_IC3.nona) 
k2=kruskal.test(IC1 ~ Grade, data = ggdata) 
A_IC3.nona=na.omit(data.frame(A_IC3,icaSetBreast_IC3$grade))
p2 <- ggplot(A_IC3.nona, aes(icaSetBreast_IC3.grade, scale(A_IC3.nona[,1])))+ geom_boxplot(aes(fill=icaSetBreast_IC3.grade),alpha=I(0.6)) + geom_jitter(colour="black",size=I(4), alpha=I(0.6))
p2=p2 + 
  
  labs(title=paste("Grade by IC1 - ICA onlycorr\n p.val=",round(k1$p.value,3)), x="Grade", y="T-cell") +
  theme_bw() +
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  #eliminates background, gridlines, and chart border
  theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_blank()
  ) +
  theme(axis.line = element_line(color = 'black'))


ggdata=na.omit(ggdata)
p2b <- ggplot( ggdata, aes(Grade, -IC1))+ geom_boxplot(aes(fill=Grade),alpha=I(0.6)) + geom_jitter(colour="black",size=I(4), alpha=I(0.6))
p2b=p2b +
  labs(title=paste("Grade by IC1 - ROMA corr\n p.val=",round(k2$p.value,3)), x="Grade", y="T-cell") +
  theme_bw() +
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  #eliminates background, gridlines, and chart border
  theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_blank()
  ) +
  theme(axis.line = element_line(color = 'black'))

grid.arrange(p2, p2b, ncol=2)

k1=kruskal.test(A_IC3.nona[,4] ~ icaSetBreast_IC3.grade, data = A_IC3.nona) 
k2=kruskal.test(IC4 ~ Grade, data = ggdata) 
A_IC3.nona=na.omit(data.frame(A_IC3,icaSetBreast_IC3$grade))
p2 <- ggplot(A_IC3.nona, aes(icaSetBreast_IC3.grade, scale(A_IC3.nona[,4])))+ geom_boxplot(aes(fill=icaSetBreast_IC3.grade),alpha=I(0.6)) + geom_jitter(colour="black",size=I(4), alpha=I(0.6))
p2=p2 + 
  
  labs(title=paste("Grade by IC4 - ICA only\n p.val=",round(k1$p.value,3)), x="Grade", y="T-cell") +
  theme_bw() +
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  #eliminates background, gridlines, and chart border
  theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_blank()
  ) +
  theme(axis.line = element_line(color = 'black'))


ggdata=na.omit(ggdata)
p2b <- ggplot( ggdata, aes(Grade, IC4))+ geom_boxplot(aes(fill=Grade),alpha=I(0.6)) + geom_jitter(colour="black",size=I(4), alpha=I(0.6))
p2b=p2b +
  labs(title=paste("Grade by IC4 - ROMA corr\n p.val=",round(k2$p.value,3)), x="Grade", y="T-cell") +
  theme_bw() +
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  #eliminates background, gridlines, and chart border
  theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_blank()
  ) +
  theme(axis.line = element_line(color = 'black'))

grid.arrange(p2, p2b, ncol=2)

p <- ggplot(ggdata, aes(icaSetBreast_IC3$age, scale(IC5)))+ geom_boxplot() + geom_jitter(colour=icaSetBreast_IC3$age)
p


summary(lm(Age ~ IC2  + IC5, data=ggdata))

plot=qplot(IC5, IC2, data=ggdata, colour=Age, size=I(10), alpha=I(0.6))
plot+theme_bw() +
  labs(title="IC2 vs. T cell vs. Age - ROMA corr\n", x="T-cell [ROMA sample projection]", y="Age [years]")+
  scale_colour_gradient(low = "blue", high="red") +
  #eliminates background, gridlines, and chart border
  theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_blank()
  ) +
  theme(axis.line = element_line(color = 'black'))

plot(ggdata$Age, ggdata$IC2)
cor.test(ggdata$Age, ggdata$IC2,method="spearman")
ct=cor.test(ggdata$Age, ggdata$IC5, method="spearman")
str(ct)

x=ggdata[,1:5]
y=ggdata$Grade
# scatterplot matrix
featurePlot(x=x, y=y, plot="ellipse")
featurePlot(x=x, y=y, plot="box")

# density plots for each attribute by class value
scales <- list(x=list(relation="free"), y=list(relation="free"))
featurePlot(x=x, y=y, plot="density", scales=scales)


#for continous varaibles
featurePlot(x = x,
            y = ggdata$Age,
            plot = "scatter",
            type = c("p", "smooth"),
            span = .5)

library(AppliedPredictiveModeling)
transparentTheme(trans = .4)

theme1 <- trellis.par.get()
theme1$plot.symbol$col = rgb(.2, .2, .2, .4)
theme1$plot.symbol$pch = 16
theme1$plot.line$col = rgb(1, 0, 0, .7)
theme1$plot.line$lwd <- 2
trellis.par.set(theme1)


#ICA AS ROMA

setwd("/Users/urszulaczerwinska/Documents/CURIE2015/Data/Biton.Zinovyev/GSE6532/ICAtest/8ics/IC3_BreastBCR_min/ICA_signif_genes")
genes_ICA_signi=read.table(file="S_significative_numerical_5.num", sep="\t", header=FALSE, row.names=NULL)[,-6]
genes_ICA_ids=read.table(file="significative_ids", sep="\t", header=FALSE, row.names=NULL)[,-6]
row.names(genes_ICA_signi)=genes_ICA_ids
s=intersect(row.names(S(icaSetBreast_IC3)),genes_ICA_ids)


genes_ICA_signi_2=read.table(file="S_significative_numerical_2_5.num", sep="\t", header=FALSE, row.names=NULL)[,-6]
row.names(genes_ICA_signi_2)=genes_ICA_ids

bf_af=data.frame(genes_ICA_signi,S(icaSetBreast_IC3)[s,])
rcorr(scale(as.matrix(bf_af)))

bf_af_2=data.frame(genes_ICA_signi_2,S(icaSetBreast_IC3)[s,])
rcorr(as.matrix(scale(bf_af_2)))


a=S(icaSetBreast_IC3)[order(-S(icaSetBreast_IC3)[,2]),][1:10,,drop=FALSE]

b=genes_ICA_signi[genes_ICA_signi[order(-genes_ICA_signi[,2]),]

b[1:10,,drop=FALSE]
