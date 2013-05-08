`break.CI`<-function(lik.bpts,pprob,bpts,my.map,index){
	
	ord<-c(1:length(bpts))
	mle.ord<-ord[lik.bpts==max(lik.bpts)]
  	denom<-lik.bpts[mle.ord]
	i<-mle.ord
	while (i < length(bpts)){
		if (lik.bpts[mle.ord]>(lik.bpts[i]+2) & pprob[i]<.01){i<-i-1;break} #considering a 2 log drop in likelihood as the CI
		i<-i+1
	}
	end<-i
	
	i<-mle.ord
	while (i > 1){
		if (lik.bpts[mle.ord]>(lik.bpts[i]+2) & pprob[i]<.01){i<-i+1;break}
		i<-i-1
	}
	start<-i

	###safeguard against multimodality
	rem<-which(pprob>=.01)
	rem<-rem[!rem %in% start:end]
	
	if (length(rem)!=0){
		start<-min(min(rem),start)
  		end<-max(max(rem),end)    
  	}
	ord = NULL
	if(index==1) {
		ord<-1:length(my.map$physposStart)
	}
	if(index==2) {
		ord<-1:length(my.map$physposEnd)
	}
	ret<-vector()
	for (i in (start:end)){
		if(index==1) {
			my.ord<-ord[my.map$physposStart==(bpts[i]-1)]
			ret[i] = cbind(as.character(my.map[my.ord,1]),as.character(my.map[my.ord+1,1]),my.map$physposStart[my.ord],my.map$physposStart[my.ord+1],round(pprob[i],digits=3))
		}
		if(index==2) {
			my.ord<-ord[my.map$physposEnd==(bpts[i]-1)]
			ret[i] = cbind(as.character(my.map[my.ord,1]),as.character(my.map[my.ord+1,1]),my.map$physposEnd[my.ord],my.map$physposEnd[my.ord+1],round(pprob[i],digits=3))
		}
	}
	return(bpts[which.max(ret)])
}
