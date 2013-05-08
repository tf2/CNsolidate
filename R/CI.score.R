`CI.score` <- function(a,b) {
	return(1-(min(a, b)/max(a, b)))
}