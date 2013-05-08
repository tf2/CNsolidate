`syn.genome` <- function(chrs) {
	rep=NULL
	data = NULL
	for(x in 1:length(chrs[,1])) {
		d = syn.data(chrs[x,1], chrs[x,2], chrs[x,3], chrs[x,4], chrs[x,5])
		data=rbind(data, d$data)
		rep=rbind(rep, d$rep)
	}
	fin = list()
	fin$data = data
	fin$rep = rep
	return(fin)
}
