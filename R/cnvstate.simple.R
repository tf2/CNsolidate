`cnvstate.simple` <- function(r) {
	s = NULL
	if(r <= -1.7) {
		s = 0
	} else if (r < 0 & r > -1.7 ) {
		s = 1
	} else if(r > 0) {
		s = 3
	}
return(s)
}