`combin.image`<- function(set=NULL) {

	file=set$files$file; odir=set$files$odir;
	n = unlist(strsplit(file, "/"))
	file1 = paste(odir, "/", substr(n[length(n)], 0, nchar(n[length(n)])-4), "_FinalReport.txt", sep="")
	file2 = paste(odir, "/", "Cmatrices/", substr(n[length(n)], 0, nchar(n[length(n)])-4), "_NumberMatrix.txt", sep="")
	file3 = paste(odir, "/", "Cmatrices/", substr(n[length(n)], 0, nchar(n[length(n)])-4), "_CombImage.jpeg", sep="")
	file4 = paste(odir, "/", "Cmatrices/", substr(n[length(n)], 0, nchar(n[length(n)])-4), "_FDR.txt", sep="")
	
	d = read.table(file1)
	m = read.table(file2, header=T)
	fdr = read.table(file4)
	m=as.matrix(m)
	m=m[d[,1]<23,]
	d=d[d[,1]<23,]
	for(x in 1:length(m[1,])) {
		m[m[,x]>1,x] = 1
	}
	s = d[,3]-d[,2]
	r=abs(d[,4])
	p = d[,5]
	s=s/max(s)
	r = r/max(r)
	p=p/max(p)
	si=apply(m,1,sum)
	m[si==0,1]=1
	l = apply(m,1,sum)
	v = vector(length=length(l))
	t = vector(length=length(l))
	for(x in 1:length(v)) {
		val = fdr[fdr[,1]==l[x],5]
		v[x] = val
		val = fdr[fdr[,1]==l[x],2]
		t[x] = val
	}
	too = t
	vo = v
	cutFdr = 5/max(v)
	v=v/max(v)
	v=v[order(l)]
	t=t/max(t)
	t=t[order(l)]
	too=too[order(l)]
	vo=vo[order(l)]
	s1= smooth.spline(1:length(l)/length(l), v, df = 5)$y
	s2= smooth.spline(1:length(l)/length(l), t, df = 5)$y
	l=l/max(l)
	m=m[order(l),]
	ll = apply(m,2,sum)
	m=m[,order(ll)]
	s=s[order(l)]
	r=r[order(l)]
	p=p[order(l)]
	l=l[order(l)]
	too=too[order(l)]
	vo=vo[order(l)]
	sss= smooth.spline(1:length(l)/length(l), l, df = 5)$y	
	st= smooth.spline(1:length(l)/length(l), r, df = length(r))$y
	sp= smooth.spline(1:length(l)/length(l), p, df = length(r))$y	
	fac = 1/(length(colnames(m))-1)
	
	jpeg(filename=file3, width=750, height=750)
	image(as.matrix(m), col=c("black", "dark grey"), axes=F)
	title(paste("Sensitivity= ", round(too[which.min(abs(s1-cutFdr))]/max(too), 4), sep=""), line=3)
	title(paste("Specificity= ", 1-(round(vo[which.min(abs(s1-cutFdr))]/100,4)), sep=""), line=2)
	title(sub=paste("Recommended number of algorithms = ", which(unique(too)==too[which.min(abs(s1-cutFdr))]), sep=""))
	axis(2, at = seq(0, 1, by = fac), labels=colnames(m))
	matplot(1:length(l)/length(l), st,add=T, type="l", col="green", lwd=2)
	matplot(1:length(l)/length(l), sp,add=T, type="l", col="yellow", lwd=2)
	#matplot(1:length(l)/length(l), sss,add=T, type="l", col="white", lwd=4)
	matplot(1:length(l)/length(l), l,add=T, type="l", col="white", lwd=4)
	matplot(1:length(l)/length(l), s1,add=T, type="l", col="red", lwd=4)
	matplot(1:length(l)/length(l), s2,add=T, type="l", col="blue", lwd=4)
	a = 1:length(l)/length(l)
	abline(v=a[which.min(abs(s1-cutFdr))], col="red", lty="dashed", lwd=4)
	dev.off()
}

