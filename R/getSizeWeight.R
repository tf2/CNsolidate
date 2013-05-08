`getSizeWeight` <- function(x) {
	c=c(4.028e-08, 5.000e-02, -1.250e-03, 2.083e-05, -2.602e-07, 2.592e-09, -2.127e-11, 1.445e-13, -7.892e-16, 3.226e-18, -8.622e-21, 1.108e-23)
	return(min(1,polyfunction(x,c)))
}
