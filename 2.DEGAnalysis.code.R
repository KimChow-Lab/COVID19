library(ggpubr)
library(ggplot2)
library(DESeq2)
library(RColorBrewer) 
library(reshape2)
setwd("D:/COVID19/GoldenHamster")

group=read.table("MetaData.txt",header=T,sep="\t",row.names=1,check.names=F,stringsAsFactors =TRUE)  
counts=read.table("FeatureCount/Hamster_counts.txt",skip=1,header=T,row.names=1) #manually changed the header

#-------------DEG analysis-------------------------------
metadata <- counts[,1:5]
countdata <- counts[,6:ncol(counts)]
countdata=countdata[,rownames(group)] #

kb <- metadata$Length / 1000
rpk <- countdata / kb 
tpm <- t(t(rpk)/colSums(rpk) * 1000000) 
tpm0=tpm[rowSums(tpm)>0,]
dim(tpm0)
#29608    40
write.table(tpm0,file="GeneExpr/GeneTPM8Symbol.txt",sep="\t",quote=F)

all(colnames(countdata)==rownames(group)) #TRUE
group$Cetegory=paste0(group$Infection,"_",group$InfectionDay,sep="")
#Control_D30, H1N1_D30, BA2_D30, COVIDWT_D30, Control_D7, H1N1_D7, BA2_D7, COVIDWT_D7
dds <- DESeqDataSetFromMatrix(countData = countdata,
                              colData = group,
                              design= ~ Cetegory)

dds <- DESeq(dds)

res <- results(dds, contrast=c("Cetegory","H1N1_D7","Control_D7"))
summary(res)
resOrdered <- res[order(res$pvalue),]
resOrdered=na.omit(resOrdered)
write.table(resOrdered,file="DEG/DEG_D7_H1N1vsCtrl.txt",quote=F,sep="\t")

sigDEG=resOrdered[resOrdered$pvalue<0.01&abs(resOrdered$log2FoldChange)>log2(1.2),]
sigDEG$Pattern=ifelse(sigDEG$log2FoldChange>0,"Up","Down")
write.table(sigDEG,file="DEG/SigDEG_D7_H1N1vsCtrl.txt",sep="\t",quote=F)


#-------------Overlap sigDEGs between different conditions-------------------------------
DEG_D7_H1N1vsCtrl=read.table("DEG/DEG_D7_H1N1vsCtrl.txt",header=T)
DEG_D30_H1N1vsCtrl=read.table("DEG/DEG_D30_H1N1vsCtrl.txt",header=T)
DEG_D7_H1N1vsCtrl$Day="D7"
DEG_D7_H1N1vsCtrl$Infection="H1N1"
DEG_D7_H1N1vsCtrl$Symbol=rownames(DEG_D7_H1N1vsCtrl)
DEG_D30_H1N1vsCtrl$Day="D30"
DEG_D30_H1N1vsCtrl$Infection="H1N1"
DEG_D30_H1N1vsCtrl$Symbol=rownames(DEG_D30_H1N1vsCtrl)

DEG_D7_BA2vsCtrl=read.table("DEG/DEG_D7_BA2vsCtrl.txt",header=T)
DEG_D30_BA2vsCtrl=read.table("DEG/DEG_D30_BA2vsCtrl.txt",header=T)
DEG_D7_BA2vsCtrl$Day="D7"
DEG_D7_BA2vsCtrl$Infection="BA2"
DEG_D7_BA2vsCtrl$Symbol=rownames(DEG_D7_BA2vsCtrl)
DEG_D30_BA2vsCtrl$Day="D30"
DEG_D30_BA2vsCtrl$Infection="BA2"
DEG_D30_BA2vsCtrl$Symbol=rownames(DEG_D30_BA2vsCtrl)

DEG_D7_COVIDWTvsCtrl=read.table("DEG/DEG_D7_COVIDWTvsCtrl.txt",header=T)
DEG_D30_COVIDWTvsCtrl=read.table("DEG/DEG_D30_COVIDWTvsCtrl.txt",header=T)
DEG_D7_COVIDWTvsCtrl$Day="D7"
DEG_D7_COVIDWTvsCtrl$Infection="COVIDWT"
DEG_D7_COVIDWTvsCtrl$Symbol=rownames(DEG_D7_COVIDWTvsCtrl)
DEG_D30_COVIDWTvsCtrl$Day="D30"
DEG_D30_COVIDWTvsCtrl$Infection="COVIDWT"
DEG_D30_COVIDWTvsCtrl$Symbol=rownames(DEG_D30_COVIDWTvsCtrl)

allDEG=rbind(DEG_D7_H1N1vsCtrl,DEG_D30_H1N1vsCtrl,DEG_D7_BA2vsCtrl,DEG_D30_BA2vsCtrl,DEG_D7_COVIDWTvsCtrl,DEG_D30_COVIDWTvsCtrl)
sigDEG=allDEG[allDEG$pvalue<0.01&abs(allDEG$log2FoldChange)>log2(1.2),]
sigDEG$Pattern=ifelse(sigDEG$log2FoldChange>0,"Up","Down")
write.table(sigDEG,file="DEG/allSigGene.txt",sep="\t",quote=F)

#Fig.1A
sigDEG$Infection=factor(sigDEG$Infection,levels=c("H1N1","BA2","COVIDWT"))
sigDEG$Day=factor(sigDEG$Day,levels=c("D7","D30"))
t=ggplot(sigDEG, aes(Infection, y=log2FoldChange,color=Pattern,size=-log10(pvalue)))+
geom_jitter(position=position_jitter(0.4),shape=1)+
ylim(-10,10)+
labs(title="",x="", y = "")+ theme_bw()+
facet_wrap(~Day)+
scale_color_manual(values=c("blue","red"))+ 
scale_size(range = c(0, 4))+
theme(axis.title.y = element_text(size=rel(1.0)),axis.text.x = element_text(size=rel(1.0)),axis.text.y = element_text(size=rel(1.0)))
pdf("Graph/sigDEG.pdf",width=8,height=4)
print(t)
dev.off()

#Fig.1F
sigDEG$Group=paste0(sigDEG$Infection,"_",sigDEG$Day,"_",sigDEG$Pattern,sep="")
sigDEGList=split(sigDEG$Symbol,sigDEG$Group)
GroupList=unique(names(sigDEGList))
edgeWeight=matrix(data = NA, nrow = length(GroupList)*(length(GroupList)-1)/2, ncol = 4)
index=0
for(i in 1:(length(GroupList)-1)){
  for (j in (i+1):length(GroupList)){
    index=index+1
    group1=sigDEGList[[i]]
    group2=sigDEGList[[j]]
    edge=length(intersect(group1,group2))/length(c(group1,group2))
    edgeWeight[index,1]=names(sigDEGList[i])
    edgeWeight[index,2]=names(sigDEGList[j])
    edgeWeight[index,3]=edge
    edgeWeight[index,4]=length(intersect(group1,group2))
  }
}
edgeWeight=data.frame(edgeWeight)
colnames(edgeWeight)=c("Node1","Node2","edgeWeight","Overlap")
edgeWeight=edgeWeight[as.numeric(edgeWeight$Overlap)>4,]
write.table(edgeWeight,file="DEG/Network/SharedDEG.txt",sep="\t",quote=F,row.names=F) #network was visualized using cytoscape



### Check the effects of normal development on H1N1 short-term effects Fig S1 ###
DEG=read.table("DEG/DEG_D7_H1N1vsCtrl.txt",header=T)
hs_data=DEG
hs_data$threshold = as.factor(ifelse(hs_data$pvalue < 0.01 & abs(hs_data$log2FoldChange) >= log2(1.2), ifelse(hs_data$log2FoldChange >= log2(1.2) ,'Up','Down'),'NoSig'))
table(hs_data$threshold)
#Down NoSig    Up 
#1385 15510   713
hs_data$ID=rownames(hs_data)
CtrlDEGByDay=read.table("DEG/DEG_Control_D30vsD7.txt",header=T,row.names=1)
CtrlSigDEGByDay=CtrlDEGByDay[CtrlDEGByDay$pvalue<0.01&abs(CtrlDEGByDay$log2FoldChange)>log2(1.2),]
CtrlSigDEGByDay$Pattern=ifelse(CtrlSigDEGByDay$log2FoldChange>0,"UpInD30","DnInD30")
hs_data$CtrlPattern =ifelse(hs_data$ID%in%rownames(CtrlSigDEGByDay),CtrlSigDEGByDay[hs_data$ID,"Pattern"],"NotShared")
table(hs_data$CtrlPattern)
#DnInD30 NotShared   UpInD30 
#240     17187       181
hs_data=hs_data[-grep("^LOC",rownames(hs_data)),]

t=ggplot(data = hs_data, aes(x = log2FoldChange, y = -log10(pvalue), size=CtrlPattern, colour=CtrlPattern, label =ID)) +
  geom_point(alpha=0.4) +
  theme_bw() + 
  scale_size_manual(values=c(3,0.5,3))+
  scale_color_manual(values=c("blue", "grey","red")) +
  #xlim(c(-5, 5)) + ylim(0,12.5)+
  geom_vline(xintercept=c(-log2(1.2),log2(1.2)),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.01),lty=4,col="black",lwd=0.8) +
  labs(x="log2 (fold change)",y="-log10 (p-value)",title="Differential Expressed Gene") +
  theme(plot.title = element_text(hjust = 0.5), legend.position="right") +  
    geom_text_repel(
    data = subset(hs_data, hs_data$pvalue < 0.01 & abs(hs_data$log2FoldChange) >= log2(1.2) & hs_data$ID%in%rownames(CtrlSigDEGByDay)),
    aes(label = ID),
    size = 3,
    max.overlaps=10,
    box.padding = unit(0.25, "lines"),
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE )
pdf("Graph/DEG_D7_H1N1.pdf",width=8,height=6)
print(t)
dev.off()





