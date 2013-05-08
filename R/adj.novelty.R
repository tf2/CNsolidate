`adj.novelty` <- function(noi=0.2, novel=0.3) {
	if(noi<0.12) {
		noi = 0.12
	}
	if(noi >0.51) {
		noi = 0.51
	}
	target=seq(0,1,by=0.001666)
	tlevel = which.min(abs(target-novel))
	data("Nmat.RData")
	p=polyfunction(noi,Nmat[tlevel,])
	return(as.numeric(p))
}
