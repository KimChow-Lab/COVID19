library(ggpubr)
library(ggplot2)
library(DESeq2)
library(RColorBrewer) 
library(reshape2)
setwd("D:/COVID19/GoldenHamster")

#-------------TPM matrix-------------------------------
group=read.table("MetaData.txt",header=T,sep="\t",row.names=1,check.names=F,stringsAsFactors =TRUE)  
counts=read.table("FeatureCount/Hamster_counts.txt",skip=1,header=T,row.names=1) #manually changed the header

kb <- metadata$Length / 1000
rpk <- countdata / kb 
tpm <- t(t(rpk)/colSums(rpk) * 1000000) 
tpm0=tpm[rowSums(tpm)>0,]
dim(tpm0)
#29608    40
write.table(tpm0,file="GeneExpr/GeneTPM8Symbol.txt",sep="\t",quote=F)


#-------------Visualization for the DEGs across conditions-------------------------------
#boxplot of target genes
GeneExpr=read.table("GeneExpr/GeneTPM8Symbol.txt",header=T,row.names=1)
dim(GeneExpr)
#29608    40
sampleInfo=read.table("MetaData.txt",header=T,sep="\t",row.names=1,check.names=F,stringsAsFactors =TRUE)  
sampleInfo$Infection=factor(sampleInfo$Infection,levels=c("Control","H1N1","BA2","COVIDWT"))
sampleInfo$InfectionDay=factor(sampleInfo$InfectionDay,levels=c("D7","D30"))
sampleInfo=sampleInfo[order(sampleInfo$Infection,sampleInfo$InfectionDay),]
GeneExpr=GeneExpr[,rownames(sampleInfo)]
all(colnames(GeneExpr)==rownames(sampleInfo))
targetGenes=c("B2m","Mnda","Rsad2","Sell")
graphDataLong=cbind(sampleInfo,t(GeneExpr[targetGenes,]))
graphDataLong=melt(graphDataLong,id=c(1:3))
graphDataLong$Infection=factor(graphDataLong$Infection,levels=c("Control","H1N1","BA2","COVIDWT"))
graphDataLong$InfectionDay=factor(graphDataLong$InfectionDay,levels=c("D7","D30"))
colnames(graphDataLong)=c(colnames(graphDataLong)[1:3],"Symbol","Expr")
graphDataLongTmp=graphDataLong
graphDataLongTmp$Group=paste0(graphDataLong$Infection,"_",graphDataLongTmp$InfectionDay,sep="")
graphDataLongTmp$Group=factor(graphDataLongTmp$Group,levels=c("Control_D7","H1N1_D7","BA2_D7","COVIDWT_D7","Control_D30","H1N1_D30","BA2_D30","COVIDWT_D30"))
graphDataLongTmp=graphDataLongTmp[graphDataLongTmp$InfectionDay=="D7",]
t=ggplot(graphDataLongTmp, aes(x=Group, y=log2(Expr+1), colour=Infection)) +
  geom_boxplot(outlier.color = "white") +
  geom_point(aes(fill=Infection), show.legend=TRUE, alpha=0.6,size=2, position = position_jitterdodge(dodge.width = 1))+
  scale_colour_brewer(palette = "Set2")+
  #scale_color_manual(values=c("#28A9A1", "#C9A77C"))+
  #scale_fill_manual(values=c("#28A9A1", "#C9A77C"))+
  theme_bw()+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.25)))+
  facet_wrap(.~Symbol,scales="free_y",ncol=2)+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x = element_text(angle=90,vjust=0.5,hjust=1))
pdf("Graph/immuneDEGD7.pdf",height=5,width=6)
print(t)
dev.off()

