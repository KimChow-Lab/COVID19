library(homologene)
library(DEGreport)
library(RColorBrewer)

#Validate the gene expression pattern in the development stage
earlyFetal=c("Early fetal","Early mid-fetal","Late mid-fetal","Late fetal","Early infancy","Late infancy","Early childhood","Middle and Late childhood")
sampleInfo=read.table("D:/Brain/GSE25219/rawData/sampleInfor.txt",header=T,row.names=1,sep="\t")
sampleInfo=sampleInfo[sampleInfo$AgePeriod%in%earlyFetal,]
sampleInfo=sampleInfo[sampleInfo$Region%in%c("MFC"),]
GeneMatrix=read.table("D:/Brain/GSE25219/rawData/GeneExprMatrix8Symbol.txt",header=T,row.names=1,sep="\t") 
sampleList=intersect(rownames(sampleInfo),colnames(GeneMatrix))
length(sampleList)
sampleInfoTmp=sampleInfo[sampleList,c("Region","Age","AgePeriod","AgeDay")]


sigGene=read.table("D:/COVID19/GoldenHamster/DEG/cellTypeExpForSigDEG.txt",header=T,sep="\t")
DownInH1N1Gene=sigGene[sigGene$Cluster=="Excitatory"&sigGene$pathway=="Down_H1N1","leadingEdge"]
DownInH1N1GeneList=strsplit(DownInH1N1Gene," ")[[1]]
DownInH1N1GeneList_Human=mouse2human(DownInH1N1GeneList)$humanGene
gene=intersect(rownames(GeneMatrix),DownInH1N1GeneList_Human) #
length(gene)
MFCExpr=GeneMatrix[gene,sampleList]
all(colnames(MFCExpr)==rownames(sampleInfoTmp))
sampleInfoTmp$AgePeriod=factor(sampleInfoTmp$AgePeriod,levels=earlyFetal)
res <- degPatterns(MFCExpr, sampleInfoTmp, time = "AgePeriod",minc = 5)
write.table(res$normalized,file="D:/COVID19/GoldenHamster/DEG/DownInH1N1_D30_NeuronDev_Pattern.txt",quote=F,sep="\t")

DownInBA2Gene=sigGene[sigGene$Cluster=="Excitatory"&sigGene$pathway=="Down_BA2","leadingEdge"]
DownInBA2GeneList=strsplit(DownInBA2Gene," ")[[1]]
DownInBA2GeneList_Human=mouse2human(DownInBA2GeneList)$humanGene
gene=intersect(rownames(GeneMatrix),DownInBA2GeneList_Human) #
length(gene)
MFCExpr=GeneMatrix[gene,sampleList]
all(colnames(MFCExpr)==rownames(sampleInfoTmp))
sampleInfoTmp$AgePeriod=factor(sampleInfoTmp$AgePeriod,levels=earlyFetal)
res <- degPatterns(MFCExpr, sampleInfoTmp, time = "AgePeriod",minc = 5)
write.table(res$normalized,file="D:/COVID19/GoldenHamster/DEG/DownInBA2_D30_NeuronDev_Pattern.txt",quote=F,sep="\t")

DownInCOVIDWTGene=sigGene[sigGene$Cluster=="Excitatory"&sigGene$pathway=="Down_COVIDWT","leadingEdge"]
DownInCOVIDWTGeneList=strsplit(DownInCOVIDWTGene," ")[[1]]
DownInCOVIDWTGeneList_Human=mouse2human(DownInCOVIDWTGeneList)$humanGene
gene=intersect(rownames(GeneMatrix),DownInCOVIDWTGeneList_Human) #
length(gene)
MFCExpr=GeneMatrix[gene,sampleList]
all(colnames(MFCExpr)==rownames(sampleInfoTmp))
sampleInfoTmp$AgePeriod=factor(sampleInfoTmp$AgePeriod,levels=earlyFetal)
res <- degPatterns(MFCExpr, sampleInfoTmp, time = "AgePeriod",minc = 5)
write.table(res$normalized,file="D:/COVID19/GoldenHamster/DEG/DownInCOVIDWT_D30_NeuronDev_Pattern.txt",quote=F,sep="\t")



###Visualization ### Fig.1J
resultTable=read.table("D:/COVID19/GoldenHamster/DEG/DownInBA2_D30_NeuronDev_Pattern.txt",header=T,sep="\t")
resultTable$AgePeriod=factor(resultTable$AgePeriod,levels=earlyFetal)
t=degPlotCluster(resultTable,time="AgePeriod",color="AgePeriod",boxes = FALSE,min_genes =3)+
geom_line(aes_string(group="genes"),alpha=0.5)+
scale_color_brewer(palette="Dark2")+theme_bw()+
theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x = element_text(size=rel(1.0),colour = "black",angle=90, vjust = 0.5, hjust=1),axis.text.y = element_text(size=rel(1.0),colour = "black"))
pdf("Graph/DownInBA2_D30_NeuronDev_Pattern_All.pdf",width=9,height=6)
print(t)
dev.off()

DnCommonPattern=resultTable[resultTable$cluster%in%c("3"),]
DnCommonPattern$AgePeriod=factor(DnCommonPattern$AgePeriod,levels=earlyFetal)
t=degPlotCluster(DnCommonPattern,time="AgePeriod",color="AgePeriod",boxes = FALSE)+
geom_line(aes_string(group="genes"),alpha=0.5)+
scale_color_brewer(palette="YlOrRd")+theme_bw()+
theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x = element_text(size=rel(1.0),colour = "black",angle=90, vjust = 0.5, hjust=1),axis.text.y = element_text(size=rel(1.0),colour = "black"))
pdf("Graph/DownInBA2_D30_NeuronDev_Pattern3.pdf",width=4,height=3)
print(t)
dev.off()

resultTable=read.table("D:/COVID19/GoldenHamster/DEG/DownInCOVIDWT_D30_NeuronDev_Pattern.txt",header=T,sep="\t")
resultTable$AgePeriod=factor(resultTable$AgePeriod,levels=earlyFetal)
t=degPlotCluster(resultTable,time="AgePeriod",color="AgePeriod",boxes = FALSE,min_genes =3)+
geom_line(aes_string(group="genes"),alpha=0.5)+
scale_color_brewer(palette="Dark2")+theme_bw()+
theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x = element_text(size=rel(1.0),colour = "black",angle=90, vjust = 0.5, hjust=1),axis.text.y = element_text(size=rel(1.0),colour = "black"))
pdf("Graph/DownInCOVIDWT_D30_NeuronDev_Pattern_All.pdf",width=7,height=6)
print(t)
dev.off()

DnCommonPattern=resultTable[resultTable$cluster%in%c("1"),]
DnCommonPattern$AgePeriod=factor(DnCommonPattern$AgePeriod,levels=earlyFetal)
t=degPlotCluster(DnCommonPattern,time="AgePeriod",color="AgePeriod",boxes = FALSE)+
geom_line(aes_string(group="genes"),alpha=0.5)+
scale_color_brewer(palette="YlOrRd")+theme_bw()+
theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x = element_text(size=rel(1.0),colour = "black",angle=90, vjust = 0.5, hjust=1),axis.text.y = element_text(size=rel(1.0),colour = "black"))
pdf("Graph/DownInCOVIDWT_D30_NeuronDev_Pattern31.pdf",width=4,height=3)
print(t)
dev.off()

