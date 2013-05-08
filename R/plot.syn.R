`plot.syn` <- function(lis) {
	
		d = lis$data
		r = lis$rep
		
		for(x in 1:length(r[,1])) {
			pad = r[x,5]
			s=which(d[,2]==r[x,2])
			ss=s-pad
			if(ss<1) { ss=1; }
			e = which(d[,3]==r[x,3])
			ee = e+pad
			if(ee>length(d[,1])) { ee=length(d[,1]); }
			plot(d[ss:ee,2], d[ss:ee,4], ylim=c(min(d[ss:ee,4]),max(d[ss:ee,4])), pch=20)
			fac = round(length(d[ss:ee,2])/5)
			if(fac%%2==0) { fac=fac+1; }
			print(fac)
			m = runmed(d[ss:ee,4],fac)
			matplot(d[ss:ee,2], m, type="l", add=T, col="purple", lwd=3)
			legend("topright", c(paste("meanRatio= ", round(r[x,4],3), sep=""),paste("numberPoints= ", r[x,5], sep=""), paste("IntervalVar= ", round(r[x,6],3), sep="")), bty="n")
			points(d[s:e,2],rep(r[x,4],r[x,5]), type="l", lwd=6, col="red")
			abline(h=0, col="blue", lwd=3)
			scan("")
		}
}
