`ws` <- function(d) { 
	vv=vector()
	for(x in 1:200) {
		vv[x] = dydLRs(d,x)
	}
return(sum(abs(diff(vv))))
}