`getRatWeight_old` <- function(x, d) {
	c=c(-2.950e-03, 7.061e-01, -1.449e+01, 2.191e+02, -1.755e+03, 8.401e+03, -2.544e+04, 4.998e+04, -6.350e+04 , 5.037e+04, -2.267e+04 , 4.423e+03)/(max(1,d*5))
	return(min(1,polyfunction(x,c)))
}
