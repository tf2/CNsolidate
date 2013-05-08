`ADM_score` <- function(d) {
	r <- .C("_SADM"
			,"d" = as.double(d[,4])
			,"dat" = as.double(d[,4])
			,"prob" = as.double(d[,5])
		 	,"size" = as.integer(length(d[,4]))
		 	,"PACKAGE" = "CNsolidate")	
	return(cbind(d[,1], d[,2], d[,3], d[,4], d[,5], r$dat, r$prob, r$d))
}