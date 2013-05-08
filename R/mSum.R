`mSum` <- function(data, win, fac) {
	st = 1
	w = win
	nlen = length(data)/win
	val = (quantile(abs(data), probs=0.68)*fac)
	cha = vector()
	for(x in 1:nlen+1) {
		seg = data[st:w] 
		cha[st:w] = c(abs(diff(cumsum(seg))),0.5)
		a <- .C("cncp", 
  				"dat" = as.double(seg), 
  				"size" = as.integer(length(seg)))
		cha[cha>val] = 1
		cha[cha!=1]=0
		st = w+1
		w=win*x	
	}
	return(cha)
}
