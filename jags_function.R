frac_model<-"model{
	for(i in 1:Ndata){
		y[i]~dnorm(mu[i], tau)
		mu[i]<-beta0[frac[i]]+beta1[frac[i]]*x[i]+z[frac[i]]
	}
	
	for(f in 1:Nfrac){
		z[f]~dnorm(mu.rand, tau.rand)
		beta0[f]~dnorm(mu0G,tau0G)
 		beta1[f]~dnorm(mu1G,tau1G)
	}
	
	tau.rand~dgamma(0.001,0.001)
	mu.rand~dnorm(0,0.0001)
	tau~dgamma(0.0001,0.0001)
	tau0~dgamma(0.0001,0.001)
 	tau1~dgamma(0.0001,0.001)
 	mu0G~dnorm(0,.0001)
	tau0G~dgamma(.0001,.0001)
	mu1G~dnorm(0,.0001)
	tau1G~dgamma(.0001,.0001)
	
	

}"
