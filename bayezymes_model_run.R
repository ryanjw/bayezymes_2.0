
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

rez<-data.frame(bayezymes_results)
names(rez)<-c("g1","g2","lm_mean","micro_mean","sm_mean","ws_mean","mean","lm_sd","micro_sd","sm_sd","ws_sd","sd","lm_hi","micro_hi","sm_hi","ws_hi","hi","lm_lo","micro_lo","sm_lo","ws_lo","lo")

write.table(rez, "bayezymes_network_results.txt",sep=" ",row.names=F)


