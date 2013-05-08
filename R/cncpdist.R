`cncpdist` <- function(data, list, dist) {
	s = vector(length=length(fea[,1])-1)
	
	for(x in 1:length(s)) {
		c = c(fea[x,3],fea[x+1,2])
		s[x] = diff(c)
	}
	
	
}
