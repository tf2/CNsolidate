`merge1` <-
function(chr, p=0.95, pp=0.2) {

	v = vector(length=length(chr[,1]))
	index = which(chr[,5]>pp)
	pin = 1;
	if (length(index) == 0) {
		index = length(chr[,1]);
	}
	
	for (x in 1:length(index)) {
		v[pin:index[x]] = mean(chr[pin:index[x],4])
		pin = index[x] + 1
	}

	v[abs(v)<quantile(abs(v), probs=p)] = 0
	my = cbind(chr, v)
		
	return(my)
}

