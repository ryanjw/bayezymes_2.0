# setwd("~/Google Drive/C-Core_no_euk_12-4-14/carbon-core-soil-paper/")
# install.packages("phyloseq")
# install.packages("plyr")
# install.packages("ggplot2")
# install.packages("reshape")
library(phyloseq)
library(plyr)
library(ggplot2)
library(reshape)
#source("http://bioconductor.org/biocLite.R")
#biocLite("phyloseq")


#

metadata = read.delim(file="./cobs_metadata.txt", row.names = 1, sep="\t", header=TRUE)
ann_data <- read.delim(sep='\t', file="./annotations_cazy.txt",header=TRUE, strip.white=TRUE, row.names=1)
data_norm <- read.delim(sep='\t', file="./cumulative-all-normreca.txt",header=TRUE, strip.white=TRUE, row.names=1)

data_norm$contig_name<-rownames(data_norm)
ann_data$contig_name<-rownames(ann_data)

merged<-merge(data_norm,ann_data,by="contig_name")

merged_sub<-subset(merged, t1 =="Bacteria" | t1=="Viruses" | t1=="Archaea")
merged_sub_fungi<-subset(merged, t2=="Fungi" )

merged_noeuk<-rbind(merged_sub,merged_sub_fungi)

# noeuk_melt<-melt(merged_noeuk[,c(1:25,27)],id="Cazy_fam" )
merged_noeuk<-merged_noeuk[,c(1,4:11,15:17,18:21,27)]

noeuk_summed<-ddply(merged_noeuk, .(Cazy_fam),summarise,
 PF_LM_H03=sum(PF_LM_H03),   PF_LM_H08=sum(PF_LM_H08),   PF_LM_H14=sum(PF_LM_H14),   PF_LM_H16=sum(PF_LM_H16),
 PF_MI_H01=sum(PF_MI_H01),   PF_MI_H06=sum(PF_MI_H06),   PF_MI_H12=sum(PF_MI_H12),   PF_MI_H13=sum(PF_MI_H13),   PF_SM_H02=sum(PF_SM_H02),
 PF_SM_H10=sum(PF_SM_H10),   PF_SM_H11=sum(PF_SM_H11),   PF_WS_H04=sum(PF_WS_H04),   PF_WS_H07=sum(PF_WS_H07),   PF_WS_H09=sum(PF_WS_H09),
PF_WS_H15=sum(PF_WS_H15) 
)

noeuk_melt<-melt(noeuk_summed, id="Cazy_fam")
noeuk_cast<-data.frame(cast(noeuk_melt, variable~Cazy_fam, value="value",fun.aggregate=sum))
SoilFrac<-c(rep("LM",4),rep("Micro",4),rep("SM",3),rep("WS",4))