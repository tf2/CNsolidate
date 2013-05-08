`lossGain`<- function(data, calls, lest=250, gest=250, lRat= -0.4, gRat = 0.3) {
	nscale = -p.scale(data, length(calls[calls[,4]<0,1]), lest, abs(lRat))
	pscale = p.scale(data, length(calls[calls[,4]>0,1]), gest, abs(gRat))
	if(nscale < -1) { nscale = -1; } 
	if(pscale >1) { pscale=1; }
	lc = calls[calls[,4]<=nscale, ];;gc = calls[calls[,4]>=pscale,]; calls=rbind(lc,gc)
	return(calls[order(calls[,1], calls[,2], calls[,3]),])
}
