`cncnps` <- function(dat, a=3, n=2, prob=0.4, dist=5000) {
	fea = NULL
	u =unique(dat[,1]) 
	co = length(dat[1,])
	for(x in 1:length(u)) {
		s=dat[dat[,1]==u[x],]
		features = NULL
		i=0; j=0; k=0; r=0;
		while(i<length(s[,1]) & r==0) {
			res <- .C("wcncp", "dat" = as.double(s[,co])
					, "prob" = as.double(prob)
					, "need" = as.integer(a)
					, "allow" = as.integer(n)
					, "index" = as.integer(0)
					, "index1" = as.integer(0)
					, "index2" = as.integer(0)
					, "size" = as.integer(length(s[,1])))
			j=res$index1
			i=res$index2
				if(j>0) {
					c=s[j:i,]
					features=rbind(features, cbind(s[j,1], s[j,2], s[i,3], median(c[,4]), length(c[,1]), mean(c[,co]))) 
				}
			k = i+1
			s=s[k:length(s[,1]),]
			if(length(s[,1])<=i) { r=1; }
		}
		
		
		
		fea=rbind(fea,features)
	}
	return(fea)
}
