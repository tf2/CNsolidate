`break.log.liks2`<-function(ldat,rdat){	
	
	lprob<-1
	rprob<-1
	if (length(ldat)>1){
		lmean<-mean(ldat)
		lsd<-sd(ldat)
		lprob<-sapply(ldat,FUN=dnorm,sd=lsd,mean=lmean)
	}	
	if (length(rdat)>1){
		rmean<-mean(rdat)
		rsd<-sd(rdat)
		rprob<-sapply(rdat,FUN=dnorm,sd=rsd,mean=rmean)
	}

return(sum(log(lprob))+sum(log(rprob)))
}
