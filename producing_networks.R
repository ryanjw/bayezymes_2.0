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

head(rez_CB)

g<-graph.edgelist(as.matrix(rez_CB[,1:2]),directed=FALSE)

# weighting the edges
E(g)$weight<-rez_CB$mean

# simplifying the network
g<-simplify(g)

# remove vertices that are not connected to the overall network
g<-delete.vertices(g,which(degree(g) < 1))


plot.igraph(g,edge.width=log1p(E(g)$weight),vertex.size=log1p(betweenness(g)),vertex.label.cex=.5)