`cSum` <- function(d, fac, s=c(50, 100, 200, 500, 1000, 2000, 5000)) {
	mat = matrix(nrow=length(d), ncol=length(s))
	for(x in 1:length(s)) {
		
		c = (ceiling(length(d)/s[x]))*s[x]-length(d)
		ss=seq(c)
		ss[1:length(ss)]=0
		m = mSum(c(d,ss),s[x], fac)
		mat[,x] = m[1:length(mat[,1])]  
	}
	return(apply(mat,1,mean))
}

  	



