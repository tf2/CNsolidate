`p.scale` <- function(x, tnum=1000, minEst = 500, fac=0.3) {
	dl = dLRs(x)
	rp = rp(x)
	Q = fac * (fac * (abs(dl-rp)*sqrt((tnum*fac)*(tnum/minEst)))+1)
	return(Q)
}