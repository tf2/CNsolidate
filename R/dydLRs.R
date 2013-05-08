`dydLRs` <-
function(x, p) {
	x=x[seq(1,length(x), p)]
return(IQR(diff(na.omit(x))) / (4 * qnorm((1 + 0.5) / 2) / sqrt(2)))
}

