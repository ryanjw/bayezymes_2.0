#  setwd("~/Google Drive/GIT/bayezymes_2.0/")
rez<-read.table("bayezymes_network_results.txt",sep=" ",check.names=F,header=T)
head(rez)
# producing master network, all of these relationships are averaged across soil fractions
rez_sub<-subset(rez, hi > 0 & lo > 0)
rez_sub_lm<-subset(rez, lm_hi > 0 & lm_lo > 0)
rez_sub_micro<-subset(rez, micro_hi > 0 & micro_lo > 0)
rez_sub_ws<-subset(rez, ws_hi > 0 & ws_lo > 0)
library(igraph)

# producing the network
g<-graph.edgelist(as.matrix(rez_sub[,1:2]),directed=FALSE)
glm<-graph.edgelist(as.matrix(rez_sub_lm[,1:2]),directed=FALSE)
gmicro<-graph.edgelist(as.matrix(rez_sub_micro[,1:2]),directed=FALSE)
gws<-graph.edgelist(as.matrix(rez_sub_ws[,1:2]),directed=FALSE)

# weighting the edges
E(g)$weight<-rez_sub$mean
E(glm)$weight<-rez_sub_lm$lm_mean
E(gmicro)$weight<-rez_sub_micro$micro_mean
E(gws)$weight<-rez_sub_ws$ws_mean

# simplifying the network
g<-simplify(g)
glm<-simplify(glm)
gmicro<-simplify(gmicro)
gws<-simplify(gws)

# remove vertices that are not connected to the overall network
g<-delete.vertices(g,which(degree(g) < 1))
glm<-delete.vertices(glm,which(degree(glm) < 1))
gmicro<-delete.vertices(gmicro,which(degree(gmicro) < 1))
gws<-delete.vertices(gws,which(degree(gws) < 1))

# What are the strongest edges (top 5)?
library(plyr)
head(arrange(rez_sub_ws, -ws_mean),5)[,1:2]

# What are the most between edges (top 5)?
df<-data.frame(get.edgelist(gws),edge.betweenness(gws,directed=F,weights=E(gws)$weight))
names(df)<-c("g1","g2","btwnness")
head(arrange(df,-btwnness),5)[,1:2]

# What are the highest degree nodes (top 5)?
df<-data.frame(degree(gws))
df<-data.frame(rownames(df),df)
names(df)<-c("node","degree")
head(arrange(df,-degree),5)[1]

# What are the highest betweenness nodes (top 5)?
df<-data.frame(betweenness(gws))
df<-data.frame(rownames(df),df)
names(df)<-c("node","btwn")
head(arrange(df,-btwn),5)[1]