`polyfunction` <- function(x,c) {
	con = c[1]; c=c[-1]
	for(i in 1:length(c)) {
		con = con + c[i]*(x^i);
	}
	return(con)
}