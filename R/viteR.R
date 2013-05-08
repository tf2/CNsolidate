`ViteR` <- function(data=NULL, states=c(-1, 0, 1), normalstate=1, ems = c(-1, 1, 0, 1, 1, 1), trans=c(0.99, 0.01, 0, 0.005, 0.99, 0.005, 0, 0.01, 0.99)) {
	if(is.null(data)) { data(test); data = test; }
	jumpy = NULL; normalstates = vector();
	for(x in 1:length(data[,1])) { normalstates[x] = normalstate; }
	u = unique(data[,1])
	for(x in 1:length(u)) {
		d = data[data[,1]==u[x],]	
		res <- .C("ViteR"
				,"data" = as.double(d[,4])
				,"states" = as.double(normalstates)
				,"emissions" = as.double(ems)
				,"transitions" = as.double(trans)
				,"dN" = as.integer(length(d[,4]))
				,"sN" = as.integer(length(states))
				,"eN" = as.integer(2)
				,"tN" = as.integer(length(states))
				,"PACKAGE" = "CNsolidate")
		jumpy = rbind(jumpy, cbind(d[,1:3], res$states))
		res = NULL;
		gc()
	}
	
	invisible(jumpy)
}

`vsegment` <- function(r=NULL, minp=3, minr=0.2) {
	ind = NULL; x = 1;
	while(x <= length(r[,1])) {
		if(r[x,4]!=2) {
			ind=c(ind,x); pin = r[x,4];
			while(r[x,4]==pin & x<=length(r[,1])) { x = x+1; }
			ind=c(ind,x); x=x-1;
		} 
		x=x+1;
	}
	p=NULL; s = seq(1,length(ind),by=2)
	for(x in 1:length(s)) {
		p=rbind(p, cbind(ind[s[x]], ind[s[x]+1]-1))
	} 
	p=p[p[,2]>p[,1],]; rr = NULL;
	for(x in 1:length(p[,1])) {
		pro = (p[x,2]-p[x,1])+1
		rr = rbind(rr,cbind(r[p[x,1],1],r[p[x,1],2], r[p[x,2],2], mean(r[p[x,1]:p[x,2],3]), pro, mean(r[p[x,1]:p[x,2],4])))
	} 
	rr=rr[abs(rr[,4])>minr & rr[,5]>minp,];
return(rr)
}