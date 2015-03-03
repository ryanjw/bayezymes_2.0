#  setwd("~/Google Drive/GIT/bayezymes_2.0/")
rez<-read.table("bayezymes_network_results.txt",sep=" ",check.names=F,header=T)

# producing master network, all of these relationships are averaged across soil fractions
rez_sub<-subset(rez, hi > 0 & lo > 0)


library(igraph)

# producing the network
g<-graph.edgelist(as.matrix(rez_sub[,1:2]),directed=FALSE)

# weighting the edges
E(g)$weight<-rez_sub$mean

# simplifying the network
g<-simplify(g)

# remove vertices that are not connected to the overall network
g<-delete.vertices(g,which(degree(g) < 1))
g

# now going to produce a network that is related to CB 3.2.1.91 (GH5,GH6,GH9,CBM2,CBM3,CBM4,CBM5,CBM10,CBM1)

dim(rez_sub)
rez_CB<-data.frame()
counter<-0
for(i in 1:dim(rez_sub)[1]){
	if((rez_sub$g1[i]=="GH5" | rez_sub$g2[i]=="GH5") | (rez_sub$g1[i]=="GH6" | rez_sub$g2[i]=="GH6") | (rez_sub$g1[i]=="GH9" | rez_sub$g2[i]=="GH9") |
	(rez_sub$g1[i]=="CBM2" | rez_sub$g2[i]=="CBM2") | (rez_sub$g1[i]=="CBM3" | rez_sub$g2[i]=="CBM3") | (rez_sub$g1[i]=="CBM4" | rez_sub$g2[i]=="CBM4" |
	(rez_sub$g1[i]=="CBM5" | rez_sub$g2[i]=="CBM5") | (rez_sub$g1[i]=="CBM10" | rez_sub$g2[i]=="CBM10") | (rez_sub$g1[i]=="CBM1" | rez_sub$g2[i]=="CBM1"))
	){rez_CB<-rbind(rez_CB,rez_sub[i,])}
	counter<-counter+1
	if(counter > 100){
		print(i/dim(rez_sub)[1])
		counter<-0
	}
}



g<-graph.edgelist(as.matrix(rez_CB[,1:2]),directed=FALSE)

# weighting the edges
E(g)$weight<-rez_CB$mean

# simplifying the network
g<-simplify(g)

# remove vertices that are not connected to the overall network
g<-delete.vertices(g,which(degree(g) < 1))


plot.igraph(g,edge.width=log1p(E(g)$weight),vertex.size=log1p(betweenness(g)),vertex.label.cex=.5)

# Incorporating metadata for plotting

ann_data <- read.delim(sep='\t', file="~/Google Drive/C-Core_no_euk_12-4-14/carbon-core-soil-paper/annotations_cazy.txt",header=TRUE, strip.white=TRUE, row.names=1)
data_norm <- read.delim(sep='\t', file="~/Google Drive/C-Core_no_euk_12-4-14/carbon-core-soil-paper/cumulative-all-normreca.txt",header=TRUE, strip.white=TRUE, row.names=1)
data_norm$contig_name<-rownames(data_norm)
ann_data$contig_name<-rownames(ann_data)

# removing non fungal euks from metadata
merged<-merge(data_norm,ann_data,by="contig_name")

merged_sub<-subset(merged, t1 =="Bacteria" | t1=="Viruses" | t1=="Archaea")
merged_sub_fungi<-subset(merged, t2=="Fungi" )

merged_noeuk<-rbind(merged_sub,merged_sub_fungi)
head(merged_noeuk)

# summing up ec and non ec enzyme families
as.vector(unique(merged_noeuk$EC))
sub1<-subset(merged_noeuk, EC == "n.d.")
sub1$EC_status<-"no_ec"
sub2<-subset(merged_noeuk, EC != "n.d.")
sub2$EC_status<-"ec"
subs<-rbind(sub1,sub2)
library(plyr)

noeuk_summed<-ddply(subs, .(Cazy_fam,Cazy_fam2,EC_status,t3),summarise,.progress="text",
 PF_LM_H03=sum(PF_LM_H03),   PF_LM_H08=sum(PF_LM_H08),   PF_LM_H14=sum(PF_LM_H14),   PF_LM_H16=sum(PF_LM_H16),
 PF_MI_H01=sum(PF_MI_H01),   PF_MI_H06=sum(PF_MI_H06),   PF_MI_H12=sum(PF_MI_H12),   PF_MI_H13=sum(PF_MI_H13),   PF_SM_H02=sum(PF_SM_H02),
 PF_SM_H10=sum(PF_SM_H10),   PF_SM_H11=sum(PF_SM_H11),   PF_WS_H04=sum(PF_WS_H04),   PF_WS_H07=sum(PF_WS_H07),   PF_WS_H09=sum(PF_WS_H09),
PF_WS_H15=sum(PF_WS_H15) 
)

head(noeuk_summed)
library(reshape)
no_euk_summed_melt<-melt(noeuk_summed, id=c("Cazy_fam","Cazy_fam2","EC_status"))

head(no_euk_summed_melt)
cazy_ec_summed<-ddply(no_euk_summed_melt, .(Cazy_fam,Cazy_fam2,EC_status),summarise, count=sum(value))

cazy_ec_cast<-data.frame(cast(cazy_ec_summed, Cazy_fam2+Cazy_fam~EC_status, value="count",fun.aggregate=sum))

cazy_ec_cast$proportion<-(cazy_ec_cast$ec/(cazy_ec_cast$ec+cazy_ec_cast$no_ec))*100
head(cazy_ec_cast)

#67001f	(0-10)
#b2182b	(11-20)
#d6604d	(21-30)
#f4a582	(31-40)
#fddbc7	(41-50)
#e0e0e0	(51-60)
#bababa	(61-70)
#878787	(71-80)
#4d4d4d	(81-90)
#1a1a1a (91-100)

#now adding a color variable based on proportion of EC
cazy_ec_cast$ec_color<-"color"
summary(cazy_ec_cast)
cols<-c("#67001f","#b2182b","#d6604d","#f4a582","#fddbc7","#e0e0e0","#bababa","#878787","#4d4d4d","#1a1a1a")

for(i in 1:dim(cazy_ec_cast)[1]){
	if(cazy_ec_cast$proportion[i] >= 0 & cazy_ec_cast$proportion[i] <= 10){cazy_ec_cast$ec_color[i]<-cols[1]}
	if(cazy_ec_cast$proportion[i] > 10 & cazy_ec_cast$proportion[i] <= 20){cazy_ec_cast$ec_color[i]<-cols[2]}
	if(cazy_ec_cast$proportion[i] > 20 & cazy_ec_cast$proportion[i] <= 30){cazy_ec_cast$ec_color[i]<-cols[3]}
	if(cazy_ec_cast$proportion[i] > 30 & cazy_ec_cast$proportion[i] <= 40){cazy_ec_cast$ec_color[i]<-cols[4]}
	if(cazy_ec_cast$proportion[i] > 40 & cazy_ec_cast$proportion[i] <= 50){cazy_ec_cast$ec_color[i]<-cols[5]}
	if(cazy_ec_cast$proportion[i] > 50 & cazy_ec_cast$proportion[i] <= 60){cazy_ec_cast$ec_color[i]<-cols[6]}
	if(cazy_ec_cast$proportion[i] > 60 & cazy_ec_cast$proportion[i] <= 70){cazy_ec_cast$ec_color[i]<-cols[7]}
	if(cazy_ec_cast$proportion[i] > 70 & cazy_ec_cast$proportion[i] <= 80){cazy_ec_cast$ec_color[i]<-cols[8]}
	if(cazy_ec_cast$proportion[i] > 80 & cazy_ec_cast$proportion[i] <= 90){cazy_ec_cast$ec_color[i]<-cols[9]}
	if(cazy_ec_cast$proportion[i] > 90 & cazy_ec_cast$proportion[i] <= 100){cazy_ec_cast$ec_color[i]<-cols[10]}
}
(unique(cazy_ec_cast$ec_color))

Vertex_order<-data.frame(V(g)$name)
names(Vertex_order)<-"Cazy_fam"
Vertex_order_color<-merge(Vertex_order, cazy_ec_cast, by="Cazy_fam")

head(Vertex_order_color)

plot.igraph(g,edge.width=log1p(E(g)$weight),vertex.size=log1p(betweenness(g)),vertex.label.cex=.5,vertex.color=Vertex_order_color$ec_color,vertex.label=NA)

# adding frame_colour for Cazy_fam2

Vertex_order_color$frame_colour<-"grey"
for(i in 1:dim(Vertex_order_color)[1]){
	if(Vertex_order_color$Cazy_fam2[i] == "CB"){Vertex_order_color$frame_colour[i]<-"#377eb8"}
	if(Vertex_order_color$Cazy_fam2[i] == "GH"){Vertex_order_color$frame_colour[i]<-"#4daf4a"}
	if(Vertex_order_color$Cazy_fam2[i] == "GT"){Vertex_order_color$frame_colour[i]<-"#984ea3"}
}

plot.igraph(g,edge.width=log1p(E(g)$weight),vertex.size=log1p(betweenness(g)),vertex.label.cex=.5,vertex.color=Vertex_order_color$ec_color,vertex.label=NA,vertex.frame.color=NA)


# making table for taxonomy

noeuk_summed<-ddply(subs, .(Cazy_fam,Cazy_fam2,EC_status,t3),summarise,.progress="text",
 PF_LM_H03=sum(PF_LM_H03),   PF_LM_H08=sum(PF_LM_H08),   PF_LM_H14=sum(PF_LM_H14),   PF_LM_H16=sum(PF_LM_H16),
 PF_MI_H01=sum(PF_MI_H01),   PF_MI_H06=sum(PF_MI_H06),   PF_MI_H12=sum(PF_MI_H12),   PF_MI_H13=sum(PF_MI_H13),   PF_SM_H02=sum(PF_SM_H02),
 PF_SM_H10=sum(PF_SM_H10),   PF_SM_H11=sum(PF_SM_H11),   PF_WS_H04=sum(PF_WS_H04),   PF_WS_H07=sum(PF_WS_H07),   PF_WS_H09=sum(PF_WS_H09),
PF_WS_H15=sum(PF_WS_H15) 
)
head(noeuk_summed)

taxa_summed_melt<-melt(noeuk_summed, id=c("Cazy_fam","Cazy_fam2","EC_status","t3"))

head(taxa_summed_melt)
taxa_ec_summed<-ddply(taxa_summed_melt, .(Cazy_fam,Cazy_fam2,EC_status,t3),summarise, count=sum(value))

taxa_ec_cast<-data.frame(cast(taxa_ec_summed, Cazy_fam2+Cazy_fam+EC_status~t3, value="count",fun.aggregate=sum))
head(taxa_ec_cast)
library(vegan)
taxa_ec_cast<-cbind(taxa_ec_cast[,1:3],decostand(taxa_ec_cast[,-c(1:3)],"total"))
head(taxa_ec_cast)

write.csv(taxa_ec_cast,"CB_network_taxonomy.csv",row.names=F)


