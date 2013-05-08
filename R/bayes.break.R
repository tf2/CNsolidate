`bayes.break`<-function(lik.bpts){
	rescale.factor<-max(lik.bpts)
	slik.bpts<-exp(lik.bpts-rescale.factor)
	denom<-sum(slik.bpts)
	pprobs<-exp(lik.bpts-rescale.factor)/denom;
return(pprobs)
}