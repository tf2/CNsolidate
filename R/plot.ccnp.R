`plot.ccnp` <- function(dat, lis) {
	u = unique(dat[,1])
	par(mfrow=c(2,1))
	for(y in 1:length(lis)) {
		d = dat[dat[,1]==u[y],4]
		cc=lis[[y]]
		st = 1
		so = 500
		s1 = st
		s2 = so
		for(x in 1:1000) {
			matplot(d[s1:s2], pch=20)
			matplot(cc[s1:s2], pch=20)
			s1 = s1+so
			s2 = s2+so
			scan("")
		}
	}
}
