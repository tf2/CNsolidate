`cncp` <- function(dat, fac=3) {
	u = unique(dat[,1])
	lis = list()
	for(x in 1:length(u)) {
		d = dat[dat[,1]==u[x],4] - median(dat[dat[,1]==u[x],4], na.rm=T)
		c1 = cSum(d,fac)
		c2 =cSum(rev(d),fac)
		cc=apply(cbind(c1,rev(c2)),1,mean)
		lis[[x]] = cc
	}
	
	return(cbind(dat, unlist(lis)))
}
