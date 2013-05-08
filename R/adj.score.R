`adj.score` <- function(score,dl) {
	cc = c(1.434e-02, 1.642e+00, 4.337e+01, -3.777e+02, 1.456e+03, -3.007e+03, 3.115e+03, -6.433e+02, -1.842e+03, 1.746e+03, -4.916e+02)
	p=polyfunction(dl,cc)
	adj.score=score-(p-(0.06-(dl/10)))
	return(adj.score)
}