#  setwd("~/Google Drive/GIT/bayezymes_2.0/")

rez<-read.table("bayezymes_network_results.txt",sep=" ",check.names=F,header=T)
head(rez)
# producing master network, all of these relationships are averaged across soil fractions
rez_sub<-subset(rez, hi > 0 & lo > 0)
rez_sub_lm<-subset(rez, lm_hi > 0 & lm_lo > 0)
rez_sub_micro<-subset(rez, micro_hi > 0 & micro_lo > 0)
rez_sub_ws<-subset(rez, ws_hi > 0 & ws_lo > 0)

library(igraph)

# producing results specific to micros
rez_CB_micro<-data.frame()
counter<-0
for(i in 1:dim(rez_sub_micro)[1]){
	if((rez_sub_micro$g1[i]=="GH5" | rez_sub_micro$g2[i]=="GH5") | (rez_sub_micro$g1[i]=="GH6" | rez_sub_micro$g2[i]=="GH6") | (rez_sub_micro$g1[i]=="GH9" | rez_sub_micro$g2[i]=="GH9") |
	(rez_sub_micro$g1[i]=="CBM2" | rez_sub_micro$g2[i]=="CBM2") | (rez_sub_micro$g1[i]=="CBM3" | rez_sub_micro$g2[i]=="CBM3") | (rez_sub_micro$g1[i]=="CBM4" | rez_sub_micro$g2[i]=="CBM4" |
	(rez_sub_micro$g1[i]=="CBM5" | rez_sub_micro$g2[i]=="CBM5") | (rez_sub_micro$g1[i]=="CBM10" | rez_sub_micro$g2[i]=="CBM10") | (rez_sub_micro$g1[i]=="CBM1" | rez_sub_micro$g2[i]=="CBM1"))
	){rez_CB_micro<-rbind(rez_CB_micro,rez_sub_micro[i,])}
	counter<-counter+1
	if(counter > 100){
		print(i/dim(rez_sub_micro)[1])
		counter<-0
	}
}

gmicro_CB<-graph.edgelist(as.matrix(rez_CB_micro[,1:2]),directed=FALSE)

# weighting the edges
E(gmicro_CB)$weight<-rez_CB_micro$mean

# simplifying the network
gmicro_CB<-simplify(gmicro_CB)

# remove vertices that are not connected to the overall network
gmicro_CB<-delete.vertices(gmicro_CB,which(degree(gmicro_CB) < 1))




rez_CB_lm<-data.frame()
counter<-0
for(i in 1:dim(rez_sub_lm)[1]){
	if((rez_sub_lm$g1[i]=="GH5" | rez_sub_lm$g2[i]=="GH5") | (rez_sub_lm$g1[i]=="GH6" | rez_sub_lm$g2[i]=="GH6") | (rez_sub_lm$g1[i]=="GH9" | rez_sub_lm$g2[i]=="GH9") |
	(rez_sub_lm$g1[i]=="CBM2" | rez_sub_lm$g2[i]=="CBM2") | (rez_sub_lm$g1[i]=="CBM3" | rez_sub_lm$g2[i]=="CBM3") | (rez_sub_lm$g1[i]=="CBM4" | rez_sub_lm$g2[i]=="CBM4" |
	(rez_sub_lm$g1[i]=="CBM5" | rez_sub_lm$g2[i]=="CBM5") | (rez_sub_lm$g1[i]=="CBM10" | rez_sub_lm$g2[i]=="CBM10") | (rez_sub_lm$g1[i]=="CBM1" | rez_sub_lm$g2[i]=="CBM1"))
	){rez_CB_lm<-rbind(rez_CB_lm,rez_sub_lm[i,])}
	counter<-counter+1
	if(counter > 100){
		print(i/dim(rez_sub_lm)[1])
		counter<-0
	}
}

glm_CB<-graph.edgelist(as.matrix(rez_CB_lm[,1:2]),directed=FALSE)

# weighting the edges
E(glm_CB)$weight<-rez_CB_lm$mean

# simplifying the network
glm_CB<-simplify(glm_CB)

# remove vertices that are not connected to the overall network
glm_CB<-delete.vertices(glm_CB,which(degree(glm_CB) < 1))


# now looking for the intersection between graphs
?graph.intersection
gint<-graph.intersection(gmicro_CB,glm_CB,byname=T,keep.all.vertices=F)
as.matrix(V(gint)$name)
plot.igraph(gint,vertex.size=log1p(betweenness(gint)),vertex.label.cex=.5)
# lots of shared stuff between the networks


# looking at specific networks with just GH48 and CBM2

rez_sub<-subset(rez, hi > 0 & lo > 0 & ((g1=="CBM2" | g2=="CBM2") | (g1=="GH48" | g2=="GH48")))
rez_sub_lm<-subset(rez, lm_hi > 0 & lm_lo > 0 & ((g1=="CBM2" | g2=="CBM2") | (g1=="GH48" | g2=="GH48")))
rez_sub_micro<-subset(rez, micro_hi > 0 & micro_lo > 0 & ((g1=="CBM2" | g2=="CBM2") | (g1=="GH48" | g2=="GH48")))
rez_sub_ws<-subset(rez, ws_hi > 0 & ws_lo > 0 & ((g1=="CBM2" | g2=="CBM2") | (g1=="GH48" | g2=="GH48")))

gmicro<-graph.edgelist(as.matrix(rez_sub_micro[,1:2]),directed=FALSE)

# weighting the edges
E(gmicro)$weight<-rez_sub_micro$mean

# simplifying the network
gmicro<-simplify(gmicro)

# remove vertices that are not connected to the overall network
gmicro<-delete.vertices(gmicro,which(degree(gmicro) < 1))

plot(gmicro)

g<-graph.edgelist(as.matrix(rez_sub[,1:2]),directed=FALSE)

# weighting the edges
E(g)$weight<-rez_sub$mean

# simplifying the network
g<-simplify(g)

# remove vertices that are not connected to the overall network
g<-delete.vertices(gmicro,which(degree(g) < 1))

plot(g)