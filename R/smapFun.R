`smapFun` <- 
function(chr1, loc) {
	
	dl = dLRs(chr1)
	rp = quantile(abs(chr1), probs=0.68, na.rm=T)
	fac = (dl + rp) / 2
	chr1 = chr1+1
	
	obs <- SMAPObservations(value = chr1, chromosome=as.character(loc[,1]), startPosition=as.numeric(loc[,2]), endPosition=as.numeric(loc[,3]), 	reporterId=as.character(1:length(chr1)))
	init.means <- c(0.4,0.7,1,1.3,1.6,3)
	for(x in 1:length(init.means)) {
		if(init.means[x]>1) {
			init.means[x] = init.means[x] + fac
		}
		if(init.means[x]<1) {
			init.means[x] = init.means[x] - fac
		}
		
	}
	init.sds <- rep(0.1,6)
	phi <- cbind(init.means, init.sds)
	hmm <- SMAPHMM(noStates=6, Phi=phi, initTrans=0.02)
	profile <- smap(hmm, obs, verbose = 0)
	probs = Q(profile)
	probs[probs==3] = 0
	probs[probs!=0] = 1
	
	return(probs)
}
