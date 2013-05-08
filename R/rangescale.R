`rangescale` <- function(Xx) {
	X=as.matrix(Xx)
	ndim <- dim(X)
	m = ndim[1]
	Xmax <- apply(X,2,max)
	Xmin <- apply(X,2,min)
	ONE <- matrix(1,m,1)
	Xscaled <- (X - (ONE %*% Xmin))/((ONE %*% Xmax) - (ONE %*% Xmin))
	return(as.vector(Xscaled))
}
