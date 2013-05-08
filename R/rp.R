`rp` <- function(x) {
	return(quantile(abs(x), probs=0.68, na.rm=T))
}
