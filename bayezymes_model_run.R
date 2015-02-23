setwd("~/Google Drive/COBS_16S_biodiversity/")
data<-read.table("phylo_and_fracs_unrar.txt",header=TRUE, check.names=FALSE, sep=" ")
setwd("~/Google Drive/postdoc_writing/C-cycling_bayezymes/16S_co-occurrence_first/")
head(data)
library(reshape)
genera<-data.frame(cast(data, Sample+Date+Block+Crop+SoilFrac~genus,value="value",fun.aggregate=sum ),check.names=FALSE)
head(genera[,1:10])
meta<-read.csv("~/Google Drive/COBS_16S_biodiversity/KBase_MGRast_Metadata_9May2013_EMB.csv",header=TRUE,check.names=FALSE)

names(meta)
meta<-meta[,c(1,6,8,12:16)]
names(meta)[1]<-"Sample"
m<-merge(genera,meta,by="Sample")

library(lme4)
library(lmerTest)
library(vegan)

m_sub<-subset(m, Crop=="CC" & Date=="12-Jul")
dim(m_sub)
names(m_sub)
m_sub[,1:10]
m_sub_data<-m_sub[,-c(1:6,327:333)]
m_sub_data2<-m_sub_data[,which(colSums(decostand(m_sub_data,"pa")) > 8)]

m_rem<-cbind(m_sub[,1:6],m_sub_data2,m_sub[,329:333])

dim(m_rem)
# m_rem[,c(1:5,145:150)]

setwd("~/Google Drive/postdoc_writing/C-cycling_bayezymes/16S_co-occurrence_first/")
write.table(m_rem, "genera_enz_for_cocur_cc_jul12.txt",row.names=FALSE,sep=" ")
summary(dataset)
dataset<-m_rem
# Using rjags library that communicates to JAGS (Just Another Gibbs Sampler)
library("rjags")


# Set up variables for the model, SoilFrac is aggregate fractions that will be random terms
Ndata<-nrow(dataset)
Nfrac<-length(unique(dataset$SoilFrac))
frac<-as.numeric(dataset$SoilFrac)

# Using the following function to implement the Bayesian hierarchical model
bayezymes<-function(var1, var2){
	x<-var1
	y<-var2
	#x<-dataset[,10]
	#y<-dataset[,11]
Ndata<-nrow(dataset)
Nfrac<-length(unique(dataset$SoilFrac))
frac<-as.numeric(dataset$SoilFrac)

# centering variables
	zx<-x-mean(x)
	zy<-y-mean(y)
#hist(zy)
# calling model that is called frac_model, this object is initialized in the script frac_model.R
# Using 10 Chains
	model_imp<-jags.model(textConnection(frac_model), data=list(
		"Ndata"=Ndata,
		"x"=zx,
		"y"=zy,
		"Nfrac"=Nfrac,
		"frac"=frac
		),n.chains=10,quiet=TRUE)

#Chains are 100000 values long and thinned every 100
	chains<-coda.samples(model=model_imp, c("beta1"),n.iter=100000,thin=100)
	#plot(chains)

	results<-matrix(nrow=0,ncol=4)
# the initial 500 values are thrown away as burn in
	for(i in 1:10){
			burned<-as.matrix(chains[[i]][501:1000,],ncol=4)
			results<-rbind(results, burned)	
	}
#head(results)
# calculating mean, median, and 95% credible interval for each variable
# also maintaining the actual values of the distribution here, standard deviation was calculated later but can be included here
	means<-c(mean(results[,1]),mean(results[,2]),mean(results[,3]),mean(results[,4]),mean(results))
	sds<-c(sd(results[,1]),sd(results[,2]),sd(results[,3]),sd(results[,4]),sd(results))
	
	high95<-c(quantile(results[,1],0.975),quantile(results[,2],0.975),quantile(results[,3],0.975),quantile(results[,4],0.975),quantile(results,0.975))
	low95<-c(quantile(results[,1],0.025),quantile(results[,2],0.025),quantile(results[,3],0.025),quantile(results[,4],0.025),quantile(results,0.025))

	result_stats<-c(means,sds,as.vector(high95),as.vector(low95))
	
	return(result_stats)

}


# results are put into the following matrix
head(dataset[,1:10])
bayezymes_results<-matrix(nrow=0,ncol=22)
for(a in 3:(dim(dataset)[2]-1)){
	
	for(b in (a+1):dim(dataset)[2]){
		#b<-6
		
		fam1<-names(dataset)[a]
		fam2<-names(dataset)[b]
		output<-c(fam1,fam2,bayezymes(dataset[,a],dataset[,b]))
		#output
		bayezymes_results<-rbind(bayezymes_results, output)
		print(a/239)
		
	}
}


## visualizing results

setwd("~/Google\ Drive/postdoc_writing/C-cycling_bayezymes/")

rez<-read.table("bayez_enz_16S_results.txt",header=TRUE,check.names=FALSE)
head(rez)
names(rez)<-c("g1","g2","lm_mean","micro_mean","sm_mean","ws_mean","mean","lm_sd","micro_sd","sm_sd","ws_sd","sd","lm_hi","micro_hi","sm_hi","ws_hi","hi","lm_lo","micro_lo","sm_lo","ws_lo","lo")

rez_sub<-subset(rez, hi > 0 & lo > 0)
library(igraph)

g<-graph.edgelist(as.matrix(rez_sub[,1:2]),directed=FALSE)
plot(g)

as.matrix((unique(rez_sub$g1)))
as.matrix((unique(rez_sub$g2)))

rez_enz<-subset(rez_sub, (g1=="AP Activity (nmol/h/g dry aggregate)" | g1=="BG Activity (nmol/h/g dry aggregate)") |
(g2=="NAG Activity (nmol/h/g dry aggregate)" | g2=="CB Activity (nmol/h/g dry aggregate)" | g2=="AP Activity (nmol/h/g dry aggregate)" | g2=="BG Activity (nmol/h/g dry aggregate)" | g2=="BX Activity (nmol/h/g dry aggregate)"))

rez_enz[,1:3]
g_enz<-graph.edgelist(as.matrix(rez_enz[,1:2]),directed=FALSE)
E(g_enz)$weight<-rez_enz$mean
plot.igraph(g_enz,edge.width=log(E(g_enz)$weight,base=2))

#looking at sub network

rez_sub<-subset(rez, sm_hi > 0 & sm_lo > 0)
rez_enz<-subset(rez_sub, (g1=="AP Activity (nmol/h/g dry aggregate)" | g1=="BG Activity (nmol/h/g dry aggregate)") |
(g2=="NAG Activity (nmol/h/g dry aggregate)" | g2=="CB Activity (nmol/h/g dry aggregate)" | g2=="AP Activity (nmol/h/g dry aggregate)" | g2=="BG Activity (nmol/h/g dry aggregate)" | g2=="BX Activity (nmol/h/g dry aggregate)"))
gmicro<-graph.edgelist(as.matrix(rez_enz[,1:2]),directed=FALSE)
E(gmicro)$weight<-rez_enz$sm_mean
plot.igraph(gmicro,edge.width=log(E(gmicro)$weight,base=2))


# now picking out species from those genera
as.matrix(unique(data$genus))

data2<-subset(data, genus=="  g__Conexibacter" |  
genus=="  g__Rhodobacter" |
genus=="  g__Anaerolinea"|
genus=="  g__Dactylosporangium"|
genus=="  g__Cohnella"|
genus=="  g__Candidatus Koribacter"|
genus=="  g__Geoderma"|
genus=="  g__OR-59"|
genus=="  g__Luteolibacter"|
genus=="  g__JG37-AG-70"|
genus=="  g__Nocardia"|
genus=="  g__Pontibacter"|
genus=="  g__Balneimonas"|
genus=="  g__Niastella"|
genus=="  g__Paenibacillus"|
genus=="  g__Ramlibacter"|
genus=="  g__Singulisphaera"|
genus=="  g__Clostridium")

species<-data.frame(cast(data2, Sample+Date+Block+Crop+SoilFrac~species,value="value",fun.aggregate=sum ),check.names=FALSE)

head(species)