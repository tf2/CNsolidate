`adj.f.t.n` <- function(noi, tlevel) {
	Tlevs = seq(1,2000,by=1);
   	data("Dmat")
	p=polyfunction(noi,Dmat[tlevel,])
	return(p)
}



`adj.f.t.n.v2` <- function(noi, tlevel=1000) {
	Tlevs = seq(1,2000,by=1);
   	data("Dmat")
	p=polyfunction(noi,Dmat[tlevel,])
	return((p/(4-noi))+(noi*0.5))
}

#`adj.f.t.n.v2` <- function(noi, tlevel=1000) {
#	Tlevs = seq(1,2000,by=1);
#   	data("Dmat")
#	p=polyfunction(noi,Dmat[tlevel,])
#	return(p/(4-noi))
#}