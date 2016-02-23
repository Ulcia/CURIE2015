setwd("/Users/urszulaczerwinska/Documents/CURIE2015/Data/DC_17.02.15")
FC.pDC=read.table("FC_T_J_pDC.txt", sep="\t", header=TRUE, row.names=1)
Activation_tumor_recognition=data.frame(read.table("ACTIVATION_TUMOR_RECOGNITION.txt", sep="\t", header=FALSE, row.names=NULL))
Immune_activation=data.frame(read.table("IMMUNE_ACTIVATION.txt", sep="\t", header=FALSE, row.names=NULL))
Immune_core=data.frame(read.table("IMMUNE_CORE.txt", sep="\t", header=FALSE, row.names=NULL))
Immune_inhibition=data.frame(read.table("IMMUNE_INHIBITION.txt", sep="\t", header=FALSE, row.names=NULL))
Immune_tumor_growth=data.frame(read.table("IMMUNE_TUMOR_GROWTH.txt", sep="\t", header=FALSE, row.names=NULL))
Immune_tumor_killing=data.frame(read.table("IMMUNE_TUMOR_KILLING.txt", sep="\t", header=FALSE, row.names=NULL))
Inhibition_tumor_recognition=data.frame(read.table("INHIBITION_TUMOR_RECOGNITION.txt", sep="\t", header=FALSE, row.names=NULL))
Migration=data.frame(read.table("MIGRATION.txt", sep="\t", header=FALSE, row.names=NULL))





FC.pDC.m_1=data.frame('Gene'=row.names(FC.pDC),'FC'=FC.pDC,Module,stringsAsFactors=FALSE)
FC.pDC.m_2=data.frame('Gene'=row.names(FC.pDC),FC.pDC,Module,stringsAsFactors=FALSE)
FC.pDC.m_3=data.frame('Gene'=row.names(FC.pDC),FC.pDC,Module,stringsAsFactors=FALSE)
FC.pDC.m_4=data.frame('Gene'=row.names(FC.pDC),FC.pDC,Module,stringsAsFactors=FALSE)
FC.pDC.m_5=data.frame('Gene'=row.names(FC.pDC),FC.pDC,Module,stringsAsFactors=FALSE)
FC.pDC.m_6=data.frame('Gene'=row.names(FC.pDC),FC.pDC,Module,stringsAsFactors=FALSE)
FC.pDC.m_7=data.frame('Gene'=row.names(FC.pDC),FC.pDC,Module,stringsAsFactors=FALSE)
FC.pDC.m_8=data.frame('Gene'=row.names(FC.pDC),FC.pDC,Module,stringsAsFactors=FALSE)





names(FC.pDC.m)[1]="Gene"

FC.pDC.m_1$Module[which(row.names(FC.pDC) %in%  Activation_tumor_recognition[,1]== TRUE)]=rep("ATR",76)

FC.pDC.m_2$Module[which(row.names(FC.pDC)  %in%  Immune_activation[,1]== TRUE)]=rep("IA",length(which(row.names(FC.pDC)  %in%  Immune_activation[,1]== TRUE)))

FC.pDC.m_3$Module[which(row.names(FC.pDC)  %in%  Immune_core[,1]== TRUE)]=rep("IC",length(which(row.names(FC.pDC)  %in%  Immune_core[,1]== TRUE)))
FC.pDC.m_4$Module[which(row.names(FC.pDC)  %in%  Immune_inhibition[,1]== TRUE)]=rep("II",length(which(row.names(FC.pDC)  %in%  Immune_inhibition[,1]== TRUE)))
FC.pDC.m_5$Module[which(row.names(FC.pDC)  %in%  Immune_tumor_growth[,1]== TRUE)]=rep("ITG",length(which(row.names(FC.pDC)  %in%  Immune_tumor_growth[,1]== TRUE)))
FC.pDC.m_6$Module[which(row.names(FC.pDC) %in%  Immune_tumor_killing[,1]== TRUE)]=rep("ITK",length(which(row.names(FC.pDC)  %in%   Immune_tumor_killing[,1]== TRUE)))
FC.pDC.m_7$Module[which(row.names(FC.pDC)  %in%  Inhibition_tumor_recognition[,1]== TRUE)]=rep("ITR",length(which(row.names(FC.pDC)  %in%   Inhibition_tumor_recognition[,1]== TRUE)))
FC.pDC.m_8$Module[which(row.names(FC.pDC)  %in%  Migration[,1]== TRUE)]=rep("M",length(which(row.names(FC.pDC) %in%  Migration[,1]== TRUE)))


FC.pDC.m_1=subset(FC.pDC.m_1,Module=="ATR")
FC.pDC.m_2=subset(FC.pDC.m_2,Module=="IA")
FC.pDC.m_3=subset(FC.pDC.m_3,Module=="IC")
FC.pDC.m_4=subset(FC.pDC.m_4,Module=="II")
FC.pDC.m_5=subset(FC.pDC.m_5,Module=="ITG")
FC.pDC.m_6=subset(FC.pDC.m_6,Module=="ITK")
FC.pDC.m_7=subset(FC.pDC.m_7,Module=="ITR")
FC.pDC.m_8=subset(FC.pDC.m_8,Module=="M")


FC.pDC.m=rbind(FC.pDC.m_1,FC.pDC.m_2,FC.pDC.m_3,FC.pDC.m_4,FC.pDC.m_5,FC.pDC.m_6,FC.pDC.m_7,FC.pDC.m_8)
boxplot(FC.pDC.m_1$Gene.1)
boxplot(FC.pDC.m$Gene.1~FC.pDC.m$Module)

library(ggplot2)
p <- ggplot(FC.pDC.m, aes(factor(Module), Gene.1))

p + geom_boxplot()

length(intersect(row.names(FC.pDC),Activation_tumor_recognition[,1]))
length(which(FC.pDC.m[,1] %in%  Activation_tumor_recognition[,1]))


DEG.FC.pDC=read.table("DEG_pDC.txt", sep="\t", header=TRUE, row.names=1)
DEG.FC.pDC.fc=DEG.FC.pDC[,1,drop=FALSE]

Module=rep("na",dim(DEG.FC.pDC.fc)[1])
DEG.FC.pDC.m_1=data.frame('Gene'=row.names(DEG.FC.pDC.fc),DEG.FC.pDC.fc,Module,stringsAsFactors=FALSE)
DEG.FC.pDC.m_2=data.frame('Gene'=row.names(DEG.FC.pDC.fc),DEG.FC.pDC.fc,Module,stringsAsFactors=FALSE)
DEG.FC.pDC.m_3=data.frame('Gene'=row.names(DEG.FC.pDC.fc),DEG.FC.pDC.fc,Module,stringsAsFactors=FALSE)
DEG.FC.pDC.m_4=data.frame('Gene'=row.names(DEG.FC.pDC.fc),DEG.FC.pDC.fc,Module,stringsAsFactors=FALSE)
DEG.FC.pDC.m_5=data.frame('Gene'=row.names(DEG.FC.pDC.fc),DEG.FC.pDC.fc,Module,stringsAsFactors=FALSE)
DEG.FC.pDC.m_6=data.frame('Gene'=row.names(DEG.FC.pDC.fc),DEG.FC.pDC.fc,Module,stringsAsFactors=FALSE)
DEG.FC.pDC.m_7=data.frame('Gene'=row.names(DEG.FC.pDC.fc),DEG.FC.pDC.fc,Module,stringsAsFactors=FALSE)
DEG.FC.pDC.m_8=data.frame('Gene'=row.names(DEG.FC.pDC.fc),DEG.FC.pDC.fc,Module,stringsAsFactors=FALSE)



DEG.FC.pDC.m_1$Module[which(row.names(DEG.FC.pDC) %in%  Activation_tumor_recognition[,1]== TRUE)]=rep("ATR",length(which(row.names(DEG.FC.pDC)  %in% Activation_tumor_recognition[,1]== TRUE)))

DEG.FC.pDC.m_2$Module[which(row.names(DEG.FC.pDC)  %in%  Immune_activation[,1]== TRUE)]=rep("IA",length(which(row.names(DEG.FC.pDC)  %in%  Immune_activation[,1]== TRUE)))

DEG.FC.pDC.m_3$Module[which(row.names(DEG.FC.pDC)  %in%  Immune_core[,1]== TRUE)]=rep("IC",length(which(row.names(DEG.FC.pDC)  %in%  Immune_core[,1]== TRUE)))
DEG.FC.pDC.m_4$Module[which(row.names(DEG.FC.pDC)  %in%  Immune_inhibition[,1]== TRUE)]=rep("II",length(which(row.names(DEG.FC.pDC)  %in%  Immune_inhibition[,1]== TRUE)))
DEG.FC.pDC.m_5$Module[which(row.names(DEG.FC.pDC)  %in%  Immune_tumor_growth[,1]== TRUE)]=rep("ITG",length(which(row.names(DEG.FC.pDC)  %in%  Immune_tumor_growth[,1]== TRUE)))
DEG.FC.pDC.m_6$Module[which(row.names(DEG.FC.pDC) %in%  Immune_tumor_killing[,1]== TRUE)]=rep("ITK",length(which(row.names(DEG.FC.pDC)  %in%   Immune_tumor_killing[,1]== TRUE)))
DEG.FC.pDC.m_7$Module[which(row.names(DEG.FC.pDC)  %in%  Inhibition_tumor_recognition[,1]== TRUE)]=rep("ITR",length(which(row.names(DEG.FC.pDC)  %in%   Inhibition_tumor_recognition[,1]== TRUE)))
DEG.FC.pDC.m_8$Module[which(row.names(DEG.FC.pDC)  %in%  Migration[,1]== TRUE)]=rep("M",length(which(row.names(DEG.FC.pDC) %in%  Migration[,1]== TRUE)))




DEG.FC.pDC.m_1=subset(DEG.FC.pDC.m_1,Module=="ATR")
DEG.FC.pDC.m_2=subset(DEG.FC.pDC.m_2,Module=="IA")
DEG.FC.pDC.m_3=subset(DEG.FC.pDC.m_3,Module=="IC")
DEG.FC.pDC.m_4=subset(DEG.FC.pDC.m_4,Module=="II")
DEG.FC.pDC.m_5=subset(DEG.FC.pDC.m_5,Module=="ITG")
DEG.FC.pDC.m_6=subset(DEG.FC.pDC.m_6,Module=="ITK")
DEG.FC.pDC.m_7=subset(DEG.FC.pDC.m_7,Module=="ITR")
DEG.FC.pDC.m_8=subset(DEG.FC.pDC.m_8,Module=="M")


DEG.FC.pDC.m=rbind(DEG.FC.pDC.m_1,DEG.FC.pDC.m_2,DEG.FC.pDC.m_3,DEG.FC.pDC.m_4,DEG.FC.pDC.m_5,DEG.FC.pDC.m_6,DEG.FC.pDC.m_7,DEG.FC.pDC.m_8)


p <- ggplot(DEG.FC.pDC.m, aes(factor(Module),  logFC ))

p + geom_boxplot()

wilcox.test(DEG.FC.pDC.m_5$logFC, DEG.FC.pDC.m_2$logFC)
mean(DEG.FC.pDC.m_1$logFC)
mean(DEG.FC.pDC.m_2$logFC)

