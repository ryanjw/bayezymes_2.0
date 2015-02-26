rez<-read.table("bayezymes_network_results.txt",sep=" ",check.names=F,header=T)

rez_sub<-subset(rez, hi > 0 & lo > 0)
library(igraph)

g<-graph.edgelist(as.matrix(rez_sub[,1:2]),directed=FALSE)
plot(g)