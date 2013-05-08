CI.plot <- function(d, m) {
	
	for(i in 1:length(m[,1])) {
	
	size = (m[i,3]- m[i,2])*2
	l1 = min(m[i,2]-size, m[i,4]-size)
	l2 = max(m[i,3]+size, m[i,5]+size)
	x=d[d[,1]==m[i,1] & d[,2] >= l1 & d[,3]<= l2,]	
	lci = CI.score(abs(m[i,2]-m[i,4]),m[i,3]- m[i,2]);
	rci = CI.score(abs(m[i,3]-m[i,5]),m[i,3]- m[i,2]);

	ma = max(abs(x[,4]))
	plot(x[,2], x[,4], pch=20, xlim=c(l1, l2),ylim=c(-ma,ma), main="Break Map")
	legend("topleft", c(paste("CI Score=", round(lci,3), sep="")), pch="")
	legend("topright", c(paste("CI Score=", round(rci,3), sep="")), pch="")
	matplot(c(m[i,2], m[i,3]), c(m[i,6], m[i,6]), col="red", lwd=3, type="l", add=T)
	abline(h=0, col="blue")
	if(m[i,2]==m[i,4]) {
		abline(v=m[i,2], col="green")
	}
	if(m[i,2]!=m[i,4]) {
		abline(v=m[i,2], col="yellow")
		abline(v=m[i,4], col="orange")
	}

	if(m[i,3]==m[i,5]) {
		abline(v=m[i,3], col="green")
	}
	if(m[i,3]!=m[i,5]) {
		abline(v=m[i,3], col="yellow")
		abline(v=m[i,5], col="orange")
	}
	scan("")
	}
}