load("~/Documents/CURIE2015/Data/Biton.Zinovyev/GSE6532/GSE6532_LUMINAL.RData")
data.tam.t=t(data.tam)
data.tam.t.n=cbind(data.frame(annot.tam[,7]),data.frame(data.tam.t))
head
head(data.tam.t.o)
dim(data.tam.t.o)

data.tam.t.o=na.omit(data.tam.t)

write.table(data.tam.t.o, file= "Breast_BCR.numerical.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
write.table(row.names(data.tam.t.o), file= "Breast_BCR.numerical_ids.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
write.table(colnames(data.tam.t.o), file= "Breast_BCR.numerical_samples.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
write.table(data.tam.t.o, file= "Breast_BCR.numerical_full.txt", quote=FALSE, sep="\t", row.names=TRUE, col.names=TRUE)

dim(annot.tam)
myrows_ann=intersect(row.names(annot.tam),row.names(data.tam.t.o))
annot.tam.o=annot.tam[myrows_ann,]
dim(annot.tam.o)


data.tam.t.n.u <- na.omit(subset(data.tam.t.n, !duplicated(data.tam.t.n[,1])))
row.names(data.tam.t.n.u)=data.tam.t.n.u[,1]

data.tam.t.n.i=data.tam.t.n[data.tam.t.n[,1] %in% datasetu.f.min[,1],]
data.tam.t.n.i <- subset(data.tam.t.n.i, !duplicated(data.tam.t.n.i[,1]))


write.table(data.tam.t.n.i[,2:278], file= "IC3_min_immune_component_numerical", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
write.table(row.names(data.tam.t.n.i[,2:278]), file= "IC3_min_immune_component_ids", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
write.table(colnames(data.tam.t.n.i[,2:278]), file= "IC3_min_immune_component_samples", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
write.table(data.tam.t.n.i[,1], file= "IC3_min_immune_component_hugo_ids", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
write.table(data.tam.t.n.i[,1:278], file= "IC3_min_immune_component.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

row.names(data.tam.t.n.i)=data.tam.t.n.i[,1]
data.tam.t.n.i.lm22=data.tam.t.n.i[sig.rows,]
setwd("/Users/urszulaczerwinska/Documents/CURIE2015/Data/Biton.Zinovyev/GSE6532/ICAtest/8ics/IC3_BreastBCR_min/CIBERSORT\ deconv/")
write.table(data.tam.t.n.i.lm22[,1:278], file= "IC3_min_immune_component_lm22genes.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

setwd("/Users/urszulaczerwinska/Documents/CURIE2015/Data/Biton.Zinovyev/GSE6532/ICAtest/8ics/IC3_BreastBCR_min/CIBERSORT\ deconv")
  
  center_apply <- function(x) {
    apply(x, 1, function(y) y - mean(y))
  }
 
  
GSE11103 <- read.delim("~/Documents/CURIE2015/Data/NewmanCIBERSORT/GSE11103_matrix_mixtures.txt", row.names=1)   
summary(GSE11103)
dim(GSE11103)
GSE11103.l=log1p(GSE11103)
GSE11103.l[1:3,1:5]

summary(GSE11103.l.c)
GSE11103.l.c[1:3,1:5]

row.mean <- apply(GSE11103.l, 1, mean)
GSE11103.l.c= sweep(GSE11103.l, 1, row.mean)


write.table(GSE11103.l.c, file= "GSE11103_numerical_all.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

write.table(row.names(GSE11103.l.c), file= "GSE11103_ids_all.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

write.table(colnames(GSE11103.l.c), file= "GSE11103_samples_all.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
write.table(GSE11103.l.c, file= "GSE11103_full_all.txt", quote=FALSE, sep="\t", row.names=TRUE, col.names=TRUE)


useMart('ensembl')
library(biomaRt)
EnsembleIDs=row.names(GSE11103.l.c)

  ensemble<-useMart("ensembl");
  hsp<-useDataset(mart=ensemble, dataset='hsapiens_gene_ensembl' );
  ids<-getBM(filters= "affy_hg_u133_plus_2",
             attributes= c("affy_hg_u133_plus_2","hgnc_id", "hgnc_symbol","description"),
             values= EnsembleIDs, mart=hsp);
GSE11103.l.c.t=cbind(row.names(GSE11103.l.c),GSE11103.l.c)
  
colnames(ids)[1]="ID" 
colnames(GSE11103.l.c.t)[1]="ID" 

total <- merge(ids,GSE11103.l.c.t,by="ID")
  

total=na.omit(total)
total.h=total[,3:16]

total.h=total.h[which(total.h$hgnc_symbol!=""),]


setwd("/Users/urszulaczerwinska/Documents/CURIE2015/Data/Abbas.2009/ICA.mix")
write.table(total.h, file= "GSE11103_hugo.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)write.table(total.h[,4:14], file= "GSE11103_numerical.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
write.table(total.h[,1], file= "GSE11103_ids.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
write.table(colnames(total.h[,4:14]), file= "GSE11103_samples.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)


GSE11103.sinatures <- read.delim("~/Documents/CURIE2015/Data/NewmanCIBERSORT/GSE11103_matrix_classes.GSE11103_matrix_pure.bm.K999.0.txt", row.names=1)
GSE11103.sinatures=GSE11103.sinatures[,1:4]

GSE11103.sinatures.s <- GSE11103.sinatures[ order(row.names(GSE11103.sinatures)), ]
head(GSE11103.sinatures.s)
GSE11103.sinatures.s.l=log1p(GSE11103.sinatures.s)

row.mean <- apply(GSE11103.sinatures.s, 1, mean)
GSE11103.sinatures.c= sweep(GSE11103.sinatures.s, 1, row.mean)

row.mean <- apply(GSE11103.sinatures.s.l, 1, mean)
GSE11103.sinatures.c.l= sweep(GSE11103.sinatures.s.l, 1, row.mean)

IC2<- read.delim("~/Documents/CURIE2015/Data/NewmanCIBERSORT/ICA2/2IC_all/S_GSE11103_numerical_all.txt_2.annot", row.names=1)   
head(IC2)
IC2=IC2[,1:2]

IC2.s <- IC2[ order(row.names(IC2)), ]
plot(IC2[,1],GSE11103.sinatures.s[,1])
dim(IC2.s)
dim(GSE11103.sinatures.s)

row.names(IC2.s)==row.names(GSE11103.sinatures.s)


intersect <- function(x, y) y[match(x, y, nomatch = 0)]

myrows=intersect(row.names(IC2.s),row.names(GSE11103.sinatures.c))
IC2.r=IC2.s[myrows,]
plot(GSE11103.sinatures.c[,1],IC2.r[,1])
plot(GSE11103.sinatures.c[,2],IC2.r[,1])
plot(GSE11103.sinatures.c[,3],IC2.r[,1])
plot(GSE11103.sinatures.c[,4],IC2.r[,1])

cor.test(GSE11103.sinatures.c[,1],IC2.r[,1])
cor.test(GSE11103.sinatures.c[,2],IC2.r[,1])
cor.test(GSE11103.sinatures.c[,3],IC2.r[,1])
cor.test(GSE11103.sinatures.c[,4],IC2.r[,1])



plot(GSE11103.sinatures.c[,1],IC2.r[,2])
plot(GSE11103.sinatures.c[,2],IC2.r[,2])
plot(GSE11103.sinatures.c[,3],IC2.r[,2])
plot(GSE11103.sinatures.c[,4],IC2.r[,2])

cor.test(GSE11103.sinatures.c[,1],IC2.r[,2])
cor.test(GSE11103.sinatures.c[,2],IC2.r[,2])
cor.test(GSE11103.sinatures.c[,3],IC2.r[,2])
cor.test(GSE11103.sinatures.c[,4],IC2.r[,2])



dim(GSE11103.sinatures.r)
dim(IC2.s)



IC3<- read.delim("~/Documents/CURIE2015/Data/NewmanCIBERSORT/ICA3/3IC_all/S_GSE11103_numerical_all.txt_3.annot", row.names=1)   
summary(IC3)
IC3=IC3[,1:3]
IC3.s <- IC3[ order(row.names(IC3)), ]

myrows=intersect(row.names(IC3.s),row.names(GSE11103.sinatures.c))


IC3.r=IC3.s[myrows,]



plot(GSE11103.sinatures.c[,1],IC3.r[,1])
plot(GSE11103.sinatures.c[,2],IC3.r[,1])
plot(GSE11103.sinatures.c[,3],IC3.r[,1])
plot(GSE11103.sinatures.c[,4],IC3.r[,1])

cor.test(GSE11103.sinatures.c[,1],IC3.r[,1],method="spearman") #0.58
cor.test(GSE11103.sinatures.c[,2],IC3.r[,1],method="spearman") #-0.53
cor.test(GSE11103.sinatures.c[,3],IC3.r[,1],method="spearman") #0.11
cor.test(GSE11103.sinatures.c[,4],IC3.r[,1],method="spearman") ##-0.028




plot(GSE11103.sinatures.c[,1],IC3.r[,2])
plot(GSE11103.sinatures.c[,2],IC3.r[,2])
plot(GSE11103.sinatures.c[,3],IC3.r[,2])
plot(GSE11103.sinatures.c[,4],IC3.r[,2])

cor.test(GSE11103.sinatures.c[,1],IC3.r[,2],method="spearman") #-0.1
cor.test(GSE11103.sinatures.c[,2],IC3.r[,2],method="spearman") ##-0.037
cor.test(GSE11103.sinatures.c[,3],IC3.r[,2],method="spearman") #0.29
cor.test(GSE11103.sinatures.c[,4],IC3.r[,2],method="spearman") #-0.37

plot(GSE11103.sinatures.c[,1],IC3.r[,3])
plot(GSE11103.sinatures.c[,2],IC3.r[,3])
plot(GSE11103.sinatures.c[,3],IC3.r[,3])
plot(GSE11103.sinatures.c[,4],IC3.r[,3])

cor.test(GSE11103.sinatures.c[,1],IC3.r[,3],method="spearman") ## -0.04
cor.test(GSE11103.sinatures.c[,2],IC3.r[,3],method="spearman") #0.23
cor.test(GSE11103.sinatures.c[,3],IC3.r[,3],method="spearman") #0.57
cor.test(GSE11103.sinatures.c[,4],IC3.r[,3],method="spearman") #-0.61




IC4<- read.delim("~/Documents/CURIE2015/Data/NewmanCIBERSORT/ICA2/4IC_all/S_GSE11103_numerical_all.txt_4.annot", row.names=1)   
summary(IC4)
IC4=IC4[,1:4]
IC4.s <- IC4[ order(row.names(IC4)), ]


myrows=intersect(row.names(IC4.s),row.names(GSE11103.sinatures.c))
IC4.r=IC4.s[myrows,]


plot(GSE11103.sinatures.c[,1],IC4.r[,1])
plot(GSE11103.sinatures.c[,2],IC4.r[,1])
plot(GSE11103.sinatures.c[,3],IC4.r[,1])

cor.test(GSE11103.sinatures.c[,1],IC4.r[,1],method="spearman") #-0.58
cor.test(GSE11103.sinatures.c[,2],IC4.r[,1],method="spearman") #0.36
cor.test(GSE11103.sinatures.c[,3],IC4.r[,1],method="spearman") #-0.1
cor.test(GSE11103.sinatures.c[,4],IC4.r[,1],method="spearman") ##0.02

plot(GSE11103.sinatures.c[,1],IC4.r[,2])
plot(GSE11103.sinatures.c[,2],IC4.r[,2])
plot(GSE11103.sinatures.c[,3],IC4.r[,2])
plot(GSE11103.sinatures.c[,4],IC4.r[,2])

cor.test(GSE11103.sinatures.c[,1],IC4.r[,2],method="spearman")##-0.027
cor.test(GSE11103.sinatures.c[,2],IC4.r[,2],method="spearman")#0.15
cor.test(GSE11103.sinatures.c[,3],IC4.r[,2],method="spearman")  #0.58
cor.test(GSE11103.sinatures.c[,4],IC4.r[,2],method="spearman") #-0.58

plot(GSE11103.sinatures.c[,1],IC4.r[,3])
plot(GSE11103.sinatures.c[,2],IC4.r[,3])
plot(GSE11103.sinatures.c[,3],IC4.r[,3])
plot(GSE11103.sinatures.c[,4],IC4.r[,3])

cor.test(GSE11103.sinatures.c[,1],IC4.r[,3],method="spearman")##-0.06
cor.test(GSE11103.sinatures.c[,2],IC4.r[,3],method="spearman")##0.19
cor.test(GSE11103.sinatures.c[,3],IC4.r[,3],method="spearman")##0.018 
cor.test(GSE11103.sinatures.c[,4],IC4.r[,3],method="spearman")# -0.35

plot(GSE11103.sinatures.c[,1],IC4.r[,4])
plot(GSE11103.sinatures.c[,2],IC4.r[,4])
plot(GSE11103.sinatures.c[,3],IC4.r[,4])
plot(GSE11103.sinatures.c[,4],IC4.r[,4])

cor.test(GSE11103.sinatures.c[,1],IC4.r[,4], method="spearman") ##0.038
cor.test(GSE11103.sinatures.c[,2],IC4.r[,4], method="spearman") #0.54
cor.test(GSE11103.sinatures.c[,3],IC4.r[,4], method="spearman") ##-0.009
cor.test(GSE11103.sinatures.c[,4],IC4.r[,4], method="spearman") #-0.5

par(mfrow=c(3,1))
plot(sort(IC3.r[,1]), pch=16 )
plot(sort(IC3.r[,2]), pch=16 )
plot(sort(IC3.r[,3]), pch=16 )


par(mfrow=c(2,2))
plot(sort(IC4.r[,1]), pch=16 )
plot(sort(IC4.r[,2]), pch=16 )
plot(sort(IC4.r[,3]), pch=16 )
plot(sort(IC4.r[,4]), pch=16 )
#-------------
A_GSE11103_numerical_all.txt_3 <- read.delim("~/Documents/CURIE2015/Data/NewmanCIBERSORT/ICA3/3IC_all/A_GSE11103_numerical_all.txt_3.num", header=FALSE)
View(A_GSE11103_numerical_all.txt_3)
A_GSE11103_numerical_all.txt_3 =A_GSE11103_numerical_all.txt_3 [,1:3]

library(gplots)
row.names(A_GSE11103_numerical_all.txt_3)=c("MixA'", "MixA''","MixA'''","MixB'","MixB''","MixB'''","MixC'","MixC''","MixC'''","MixD'","MixD''","MixD'''")
colnames(A_GSE11103_numerical_all.txt_3)= c("ICA","IC2","IC3")

heatmap.2(as.matrix(A_GSE11103_numerical_all.txt_3), dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none')


A_GSE11103_numerical_all.txt_2 <- read.delim("~/Documents/CURIE2015/Data/NewmanCIBERSORT/ICA3/2IC_all/A_GSE11103_numerical_all.txt_2.num", header=FALSE)
A_GSE11103_numerical_all.txt_2 =A_GSE11103_numerical_all.txt_2 [,1:2]
row.names(A_GSE11103_numerical_all.txt_2)=c("MixA'", "MixA''","MixA'''","MixB'","MixB''","MixB'''","MixC'","MixC''","MixC'''","MixD'","MixD''","MixD'''")
colnames(A_GSE11103_numerical_all.txt_2)= c("ICA","IC2")
heatmap.2(as.matrix(A_GSE11103_numerical_all.txt_2), dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none')

A_GSE11103_numerical_all.txt_4 <- read.delim("~/Documents/CURIE2015/Data/NewmanCIBERSORT/ICA3/4IC_all/A_GSE11103_numerical_all.txt_4.num", header=FALSE)
A_GSE11103_numerical_all.txt_4 =A_GSE11103_numerical_all.txt_4 [,1:4]
row.names(A_GSE11103_numerical_all.txt_4)=c("MixA'", "MixA''","MixA'''","MixB'","MixB''","MixB'''","MixC'","MixC''","MixC'''","MixD'","MixD''","MixD'''")
colnames(A_GSE11103_numerical_all.txt_4)= c("ICA","IC2","IC3","IC4")
heatmap.2(as.matrix(A_GSE11103_numerical_all.txt_4), dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none')


#-------------------------------

esGolub$Sample <- NULL
res <- nmf(esGolub, 3)
res
fit(res)
V.hat <- fitted(res)
dim(V.hat)
summary(res)
esGolub$Cell

GSE11103.short=GSE11103[1:200, ]
GSE11103.nmf= nmf(GSE11103.short,3, .options="t")
summary(GSE11103.nmf)

layout(1)
# basis components
basismap(GSE11103.nmf, subsetRow = TRUE)
# mixture coefficients
coefmap(GSE11103.nmf)

GSE11103.nmf.4= nmf(GSE11103.short,4, .options="t")


coefmap(GSE11103.nmf.4)
summary(GSE11103.nmf.4)
coef(GSE11103.nmf)
heatmap(basis(GSE11103.nmf))

GSE11103.nmf.5= nmf(GSE11103.short,5, .options="t")
coefmap(GSE11103.nmf.5)

#--------
GSE11103.med=GSE11103[1:4000, ]
GSE11103.nmf.3= nmf(GSE11103,3, .options="t")
coefmap(GSE11103.nmf.3)

GSE11103.nmf.a.4= nmf(GSE11103,4, .options="t")
coefmap(GSE11103.nmf.a.4)
basismap(GSE11103.nmf.a.4)


emix.nmf.4= nmf(exprs(emix),4, .options="t")
profplot(coef(emix),coef(emix.nmf.4))
coef.emix.nmf=coef(emix.nmf.4)
profplot(coef(emix),coef.emix.nmf[c(2,3,4,1),])

cor.test(coef(emix)[3,],coef.emix.nmf[4,])
aheatmap(basis(emix.nmf.4))
#----
emix.nmf.4.mix= nmf(exprs(emix)[,13:24],4, .options="t")
coef.emix.nmf.mix=coef(emix.nmf.4.mix)
coef.emix=coef(emix)[,13:24]
profplot(coef.emix,coef.emix.nmf.mix[c(2,4,3,1),])

cor.test(coef.emix[4,],coef.emix.nmf.mix[1,],method="spearman")

#--------------------
#PCA
library(FactoMineR)
res.pca = PCA(GSE11103.l.c, scale.unit=TRUE, ncp=5, graph=T)
res.pca = PCA(t(GSE11103.l.c), scale.unit=TRUE, ncp=5, graph=T)


library(ggplot2)
library(grid)

PCbiplot2 <- function(res.pca, x="Dim.1", y="Dim.2") {
  if(!require(ggplot2)) install.packages("ggplot2")
  # res.pca being a PCA object
  data <- data.frame(obsnames=row.names(res.pca$ind$coord), res.pca$ind$coord)
  plot <- ggplot(data, aes_string(x=x, y=y)) + geom_text(alpha=.4, size=3,     aes(label=obsnames))
  plot <- plot + geom_hline(aes(0), size=.2) + geom_vline(aes(0), size=.2)
  datapc <- data.frame(varnames=rownames(res.pca$var$coord), res.pca$var$coord)
  mult <- min(
    (max(data[,y]) - min(data[,y])/(max(datapc[,y])-min(datapc[,y]))),
    (max(data[,x]) - min(data[,x])/(max(datapc[,x])-min(datapc[,x])))
  )
  datapc <- transform(datapc,
                      v1 = .7 * mult * (get(x)),
                      v2 = .7 * mult * (get(y))
  )
  plot <- plot + coord_equal() + geom_text(data=datapc, aes(x=v1, y=v2,     label=varnames), size = 5, vjust=1, color="red")
  plot <- plot + geom_segment(data=datapc, aes(x=0, y=0, xend=v1, yend=v2),     arrow=arrow(length=unit(0.2,"cm")), alpha=0.75, color="red")
  plot
}
PCbiplot2(res.pca)

#---------------------
setwd("/Users/urszulaczerwinska/Documents/CURIE2015/Data/NewmanCIBERSORT/ICA_MixA")


write.table(GSE11103.l.c[,1:3], file= "GSE11103_numerical_all.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

write.table(row.names(GSE11103.l.c[,1:3]), file= "GSE11103_ids_all.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

write.table(colnames(GSE11103.l.c[,1:3]), file= "GSE11103_samples_all.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
write.table(GSE11103.l.c[,1:3], file= "GSE11103_full_all.txt", quote=FALSE, sep="\t", row.names=TRUE, col.names=TRUE)

IC2_MixA<- read.delim("~/Documents/CURIE2015/Data/NewmanCIBERSORT/ICA_MixA/IC2/S_GSE11103_numerical_all.txt_2.annot", row.names=1)   
head(IC2_MixA)
colnames(IC2_MixA)=c("IC1", "IC2")
IC2_MixA=IC2_MixA[,1:2]

IC2.s <- IC2_MixA[ order(row.names(IC2_MixA)), ]
plot(IC2[,1],GSE11103.sinatures.s[,1])
dim(IC2.s)
dim(GSE11103.sinatures.s[,1:2])

row.names(IC2.s)==row.names(GSE11103.sinatures.s[,1:2])

#
intersect <- function(x, y) y[match(x, y, nomatch = 0)]

myrows=intersect(row.names(IC2.s),row.names(GSE11103.sinatures.c))
IC2.r=IC2.s[myrows,]
plot(GSE11103.sinatures.c[,1],IC2.r[,1])
plot(GSE11103.sinatures.c[,2],IC2.r[,1])
plot(GSE11103.sinatures.c[,3],IC2.r[,1])
plot(GSE11103.sinatures.c[,4],IC2.r[,1])

cor.test(GSE11103.sinatures.c[,1],IC2.r[,1],method="spearman") ##0.05
cor.test(GSE11103.sinatures.c[,2],IC2.r[,1],method="spearman") #0.172
cor.test(GSE11103.sinatures.c[,3],IC2.r[,1],method="spearman") ##0.08
cor.test(GSE11103.sinatures.c[,4],IC2.r[,1],method="spearman") ##-0.03



plot(GSE11103.sinatures.c[,1],IC2.r[,2])
plot(GSE11103.sinatures.c[,2],IC2.r[,2])
plot(GSE11103.sinatures.c[,3],IC2.r[,2])
plot(GSE11103.sinatures.c[,4],IC2.r[,2])

cor.test(GSE11103.sinatures.c[,1],IC2.r[,2],method="spearman") #0.54
cor.test(GSE11103.sinatures.c[,2],IC2.r[,2],method="spearman") #-0.56
cor.test(GSE11103.sinatures.c[,3],IC2.r[,2],method="spearman") #-0.19
cor.test(GSE11103.sinatures.c[,4],IC2.r[,2],method="spearman") #0.26


IC3_MixA<- read.delim("~/Documents/CURIE2015/Data/NewmanCIBERSORT/ICA_MixA/IC3/S_GSE11103_numerical_all.txt_3.annot", row.names=1)   
summary(IC3_MixA)
IC3_MixA=IC3_MixA[,1:3]
colnames(IC3_MixA)=c("IC1","IC2", "IC3")
IC3.s <- IC3_MixA[ order(row.names(IC3_MixA)), ]

myrows=intersect(row.names(IC3.s),row.names(GSE11103.sinatures.c))
IC3.r=IC3.s[myrows,]
plot(GSE11103.sinatures.c[,1],IC3.r[,1])
plot(GSE11103.sinatures.c[,2],IC3.r[,1])
plot(GSE11103.sinatures.c[,3],IC3.r[,1])
plot(GSE11103.sinatures.c[,4],IC3.r[,1])

cor.test(GSE11103.sinatures.c[,1],IC3.r[,1],method="spearman") #0.547
cor.test(GSE11103.sinatures.c[,2],IC3.r[,1],method="spearman") #-0.5649
cor.test(GSE11103.sinatures.c[,3],IC3.r[,1],method="spearman") #-0.19
cor.test(GSE11103.sinatures.c[,4],IC3.r[,1],method="spearman") #0.26




plot(GSE11103.sinatures.c[,1],IC3.r[,2])
plot(GSE11103.sinatures.c[,2],IC3.r[,2])
plot(GSE11103.sinatures.c[,3],IC3.r[,2])
plot(GSE11103.sinatures.c[,4],IC3.r[,2])

cor.test(GSE11103.sinatures.c[,1],IC3.r[,2],method="spearman") #-0.14
cor.test(GSE11103.sinatures.c[,2],IC3.r[,2],method="spearman") #-0.11
cor.test(GSE11103.sinatures.c[,3],IC3.r[,2],method="spearman") ##0.065
cor.test(GSE11103.sinatures.c[,4],IC3.r[,2],method="spearman") ##0.03

plot(GSE11103.sinatures.c[,1],IC3.r[,3])
plot(GSE11103.sinatures.c[,2],IC3.r[,3])
plot(GSE11103.sinatures.c[,3],IC3.r[,3])
plot(GSE11103.sinatures.c[,4],IC3.r[,3])

cor.test(GSE11103.sinatures.c[,1],IC3.r[,3],method="spearman") ## 0.10
cor.test(GSE11103.sinatures.c[,2],IC3.r[,3],method="spearman") ##-0.06
cor.test(GSE11103.sinatures.c[,3],IC3.r[,3],method="spearman") ##-0.022
cor.test(GSE11103.sinatures.c[,4],IC3.r[,3],method="spearman") #0.026


#-----------
#only signature genes
myrows=intersect(row.names(GSE11103.l.c),row.names(GSE11103.sinatures.c))
GSE11103.l.c.sig=GSE11103.l.c[myrows,]

#PCA
library(FactoMineR)
res.pca = PCA(GSE11103.l.c.sig, scale.unit=TRUE, ncp=5, graph=T)

res.pca = PCA(GSE11103.sinatures.c, scale.unit=TRUE, ncp=5, graph=T)
dim(GSE11103.sinatures.c)

res.pca$ind$coord
round(res.pca$ind$contrib,2)
res.pca
barplot(res.pca$eig[,1],main="Eigenvalues",names.arg=1:nrow(res.pca$eig))
res.pca$ind$coord
res.pca$ind$cos2
res.pca$ind$contrib
dimdesc(res.pca,axes = c(1,2))
print(res.pca)


mega=cbind(GSE11103.l.c.sig,GSE11103.sinatures.c.l)

res.pca <- PCA(mega, quanti.sup =13:16 , graph=FALSE)


res.pca$quanti.sup[,1:2]



#------------CellMix
if( !require(BiocInstaller) ){
  # enable Bioconductor repositories
  # -> add Bioc-software
  setRepositories() 
  
  install.packages('BiocInstaller')
  library(BiocInstaller)
}
biocLite("GEOquery")
install.packages("devtools")
library(devtools)
install_github("fawda123/ggord")
library(ggord)
library(CellMix)
library(GEOquery)
emix <- ExpressionMix("GSE11058")

coef(emix)[, 13:24]
head(basis(emix), 3)
dim(head(exprs(emix)))
kl <- ged(emix, basis(emix), "lsfit", verbose=TRUE)
round(coef(kl),2)
identical(basis(kl), basis(emix))

kl <- ged(emix, basis(emix), "lsfit", fit="nnls",verbose=TRUE)
identical(basis(kl), basis(emix))

profplot(coef(emix),coef(kl))

kl2 <- ged(as.matrix(GSE11103), basis(emix), "lsfit", verbose=TRUE)
profplot(coef(emix)[, 13:24],coef(kl2))


fit=lsfit(basis(emix), exprs(emix))
exprs(emix)
dim(as.matrix(GSE11103))
dim(basis(emix))
dim(exprs(emix))

round(fit$coef,2)


profplot(coef(emix),fit$coef)

qprog=ged(emix, basis(emix), "qprog")
coef(qprog)
profplot(coef(emix),coef(qprog))

deconf=ged(exprs(emix), 4, "deconf")
round(coef(deconf),2)

cor.test(coef(emix)[3,],coef(deconf)[1,])
coef.deconf=rbind(coef(deconf)[2,],coef(deconf)[1,],coef(deconf)[3,],coef(deconf)[4,])
profplot(coef(emix),coef.deconf)
heatmap(basis(deconf))
heatmap(as.matrix(GSE11103.sinatures.s))
heatmap(as.matrix(basis(emix)))

#-----Basis matrix from pdf in the paper

Abbas.basisMatrix <- read.csv("~/Documents/CURIE2015/Data/Abbas.2009/Abbas.basisMatrix.csv", header=FALSE, row.names=1, sep=";")
summary(Abbas.basisMatrix)
dim(Abbas.basisMatrix)
Abbas.basisMatrix.num=Abbas.basisMatrix[,3:19]
dim(Abbas.basisMatrix.num)
summary(Abbas.basisMatrix.num)
heatmap(as.matrix(Abbas.basisMatrix.num))

#fail - data for leukocytes

#---basis for Abbas generated by cibersort

GSE11103_matrix_classes <- read.table("~/Documents/CURIE2015/Data/NewmanCIBERSORT/GSE11103_matrix_classes.GSE11103_matrix_pure.bm.K999.0 (1).txt", header=TRUE, dec=".", row.names=1)

dim(GSE11103_matrix_classes)
GSE11103_matrix_classes=GSE11103_matrix_classes[,1:4]
summary(GSE11103_matrix_classes)
heatmap(as.matrix(-log2(GSE11103_matrix_classes)))


data(Abbas)
Abbas
exprs(Abbas)
dim(exprs(Abbas))

heatmap(log(exprs(Abbas)))


#---------HEATMAP--------------------------------
# random data that follow an 3-rank NMF model (with quite some noise: sd=2)
X <- syntheticNMF(100, 3, 20, noise = 2)
head(X)
n
# row annotations and covariates
n <- nrow(X)
d <- rnorm(n)
e <- unlist(mapply(rep, c("X", "Y", "Z"), 10))
e <- c(e, rep(NA, n - length(e)))
rdata <- data.frame(Var = d, Type = e)
# column annotations and covariates
p <- ncol(X)
a <- sample(c("alpha", "beta", "gamma"), p, replace = TRUE)

# define covariates: true groups and some numeric variable
c <- rnorm(p)
# gather them in a data.frame
covariates <- data.frame(a, X$pData, c)

par(mfrow = c(1, 2))
aheatmap(X, annCol = covariates, annRow = X$fData)
aheatmap(X)

res <- nmf(X, 3, nrun = 10)
res
opar <- par(mfrow = c(1, 2))
# coefmap from multiple run fit: includes a consensus track
coefmap(res)
# coefmap of a single run fit: no consensus track
coefmap(minfit(res))
par(opar)
coefmap(res, Colv = 'euclidean', main = "Metagene contributions in each sample", labCol = NULL, annRow = list(Metagene=':basis'), annCol = list(':basis', Class=a, Index=c), annColors = list(Metagene='Set2'), info = TRUE)


#---------------------------------
#finding_marker_genes
emix.pure=pureSamples(emix) 
emix.markers=extractMarkers(emix.pure, emix.pure$Type)
emix.markers=extractMarkers(emix.pure, emix.pure$Type)

summary(emix.markers)
hist(emix.markers,breaks=20)
breakdown(emix.markers)
barplot(emix.markers)
show(emix.markers)
emix.markers.c=emix.markers<=10^-5
summary(emix.markers)
hist(emix.markers.c)

head(geneIds(emix.markers)$"IM-9")

sel=screeplot(emix.markers, basis(emix), range=1:1500)
summary(sel)
marknames(sel)
cbind(geneValues(emix.markers)$Jurkat,geneValues(emix.markers)$"IM-9")

#heatmap(as.matrix(cbind(geneValues(emix.markers)$Jurkat,geneValues(emix.markers)$"IM-9")))

aheatmap(emix, annRow = sel$Type)

geneIds(sel)$'IM-9'
rows=c(rep("IM-9",length(geneIds(sel)$'IM-9')),rep("Jurkat",length(geneIds(sel)$'Jurkat')),rep("Raji",length(geneIds(sel)$'Raji')),rep("THP-1",length(geneIds(sel)$'THP-1')))
markergenes.matrix=exprs(emix)[c(geneIds(sel)$'IM-9',geneIds(sel)$'Jurkat',geneIds(sel)$'Raji',geneIds(sel)$'THP-1'),]
aheatmap(log(markergenes.matrix[,1:12]), annRow = as.factor(rows), Rowv=NA,annColors="rainbow")

sel.m =MarkerList(sel)

markermap(sel.m,log(markergenes.matrix[,1:12]), subsetRow = TRUE)

markermap(sel.m, emix.pure[2000:2500,], view='single')


emix.nmf.4.mark= nmf(markergenes.matrix[,13:24],4, .options="t")
profplot(coef(emix)[,13:24],coef(emix.nmf.4.mark)[c(1,4,3,2),])
coef.emix.nmf=coef(emix.nmf.4.mark)
profplot(coef(emix),coef.emix.nmf[c(4,2,1,3),])

cor.test(coef(emix)[3,13:24],coef.emix.nmf[2,])

basismap(emix.nmf.4.mark,annRow = as.factor(rows),Rowv=NA)
basismap(emix.nmf.4.mark,annRow = as.factor(rows),Rowv=NA)

basismap(emix.nmf.4.mark, subsetRow=FALSE)
aheatmap(t(scale(t(basis(emix.nmf.4.mark)))), annRow = as.factor(rows), Rowv=NA,annColors="rainbow")


basismarkermap(sel.m, emix.nmf.4.mark, annRow = as.factor(rows),Rowv=NA)

markermap(sel, emix.nmf.4.mark.r, sel$Type , subsetRow = TRUE)

basis(emix.nmf.4.mark)[,1]
myrows.m=intersect(row.names(emix.nmf.4.mark),row.names(exprs(emix)[flatten(geneIds(sel)),1:12]))
emix.nmf.4.mark.r=basis(emix.nmf.4.mark)[myrows.m,]

emix.nmf.4.mark.r=t(scale(t(basis(emix.nmf.4.mark)[myrows.m,])))
heatmap(log2(markergenes.matrix[,13:24]))
library(gplots)
heatmap.2((log(basis(emix.nmf.4.mark))))

plot(emix.nmf.4.mark.r[,3],exprs(emix)[flatten(geneIds(sel)),1:12][,11])


cor.test(emix.nmf.4.mark.r[,3],exprs(emix)[flatten(geneIds(sel)),1:12][,11],method='spearman')


head(basis(emix))
markerScoreMethod()

write.table(markergenes.matrix , file= "markergenes.matrix_numerical", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
write.table(row.names(markergenes.matrix), file= "markergenes.matrix_ids.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
write.table(colnames(markergenes.matrix), file= "markergenes.matrix_samples.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)


#--------------------marker genes ICA-------------------------------------------

mix.markers <- read.delim("~/Documents/CURIE2015/Data/Abbas.2009/ICA/S_markergenes.matrix_numerical.annot", row.names=1)   
summary(mix.markers)
head(mix.markers)
colnames(mix.markers)=c("IC1","IC2", "IC3", "IC4")
mix.markers=mix.markers[,1:4]
mix.markers.s=t(scale(t(mix.markers)))
markermap(sel, mix.markers.s, sel$Type , subsetRow = TRUE)

install.packages("Hmisc")
library(Hmisc)

emix.markers=basis(emix)[row.names(mix.markers.s),]
markers.w=cbind(mix.markers.s,emix.markers)
markers.cor=rcorr(markers.w, type="spearman")
markers.cor$P<0.05

# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

mcor=flattenCorrMatrix(markers.cor$r, markers.cor$P)

mcor.i=mcor[c(7,8,9,10,11,12,13,14,16,17,18,19, 22,23,24,25),]
mcor.i$p=round(mcor.i$p,2)
  
install.packages("corrplot")
library(corrplot)
markers.i=markers.cor$r[1:4,5:8]
corrplot(markers.i, type="upper", order="original", tl.col="black", tl.srt=45, col=colmy(300))


col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061"))
vec_col=c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061")

colmy <- colorRampPalette(rev(vec_col))

mix.samples <- read.delim("~/Documents/CURIE2015/Data/Abbas.2009/ICA/A_markergenes.matrix_numerical_4.num", row.names=NULL, header=FALSE)   
mix.samples=mix.samples [,1:4]
colnames(mix.samples)=c("IC1","IC2", "IC3", "IC4")

colnames(emix)=row.names(mix.samples)

mix.samples.s=scale(mix.samples)
aheatmap(mix.samples.s, Rowv=NULL)
profplot(coef(emix),t(mix.samples.s))

#---mixes only
write.table(markergenes.matrix[,13:24], file= "markergenes.matrix.mix_numerical", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
write.table(row.names(markergenes.matrix[,13:24]), file= "markergenes.matrix.mix_ids.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
write.table(colnames(markergenes.matrix[,13:24]), file= "markergenes.matrix.mix_samples.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)


mix.markers.mix <- read.delim("~/Documents/CURIE2015/Data/Abbas.2009/ICA.mix/S_markergenes.matrix.mix_numerical.annot", row.names=1)   
summary(mix.markers.mix)
head(mix.markers)
colnames(mix.markers.mix)=c("IC1","IC2", "IC3", "IC4")
mix.markers.mix=mix.markers.mix[,1:4]
mix.markers.mix.s=t(scale(t(mix.markers.mix)))
markermap(sel, mix.markers.mix.s, annRow=rows , subsetRow = TRUE)


emix.markers=basis(emix)[row.names(mix.markers.mix.s),]
markers.mix.w=cbind(mix.markers.mix.s,emix.markers)
markers.mix.cor=rcorr(markers.mix.w, type="spearman")

mcor.mix=flattenCorrMatrix(markers.mix.cor$r, markers.mix.cor$P)

mcor.mix.i=mcor.mix[c(7,8,9,10,11,12,13,14,16,17,18,19, 22,23,24,25),]
mcor.mix.i$p=round(mcor.mix.i$p,2)
aheatmap(mix.markers.mix.s, Rowv=NULL, annRow=rows)

markers.mix.i=markers.mix.cor$r[1:4,5:8]

corrplot(markers.mix.i, type="full", order="original", tl.col="black", tl.srt=45, col=colmy(300))

mix.samples.mix <- read.delim("~/Documents/CURIE2015/Data/Abbas.2009/ICA.mix/A_markergenes.matrix.mix_numerical_4.num", row.names=NULL, header=FALSE)   
mix.samples.mix=mix.samples.mix [,1:4]
colnames(mix.samples.mix)=c("IC1","IC2", "IC3", "IC4")

row.names(mix.samples.mix)=colnames(emix)[13:24]

mix.samples.mix.s=scale(mix.samples.mix)
aheatmap(mix.samples.mix.s, Rowv=NULL)
profplot(coef(emix)[,13:24],t(mix.samples.mix.s)[c("IC2","IC4","IC1","IC3"),])

dev.off()

setwd("/Users/urszulaczerwinska/Documents/CURIE2015/Data/Biton.Zinovyev/GSE6532/ICAtest/8ics/IC3_BreastBCR_min/Immune_sig_gmt")
Bindea=read.table("/Users/urszulaczerwinska/Documents/CURIE2015/Data/Biton.Zinovyev/GSE6532/ICAtest/8ics/IC3_BreastBCR_min/Immune_sig_gmt/Bindea_markers_2.txt", sep="\t", headers=TRUE)
for (i in 1:24){



write.table(paste(t(unique(Bindea[,i],sep="\t"))), file=paste("Bindea_unique_",i,".gmt", sep=""), quote=FALSE, sep="\t", eol = "\t",col.names=FALSE, row.names=FALSE)              

}            
               