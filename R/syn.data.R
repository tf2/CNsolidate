`syn.data` <- function(chr, l, m, sd, pro) {
		
		vars = seq(0.1,1,by=0.05) 
		rats = c(0.5,1,1.5,-0.5,-1,-2,-3) 
		lens =c(seq(3,50,by=1),100,200,500,1000)
		spacing = c(seq(1,200,by=2),500,seq(1000,10000,by=1000))
		
		lis = list()
		d = rnorm(l,m,sd)
		many = (l/100)*pro
		s = seq(1,l, by=l/many)
		rep = NULL
		for(y in 1:length(s)) {
			r = c(runif(1,1,length(lens)),runif(1,1,length(rats)), runif(1,1,length(vars)))
			r1 = rnorm(lens[r[1]],rats[r[2]],vars[r[3]])
			if(rats[r[2]]>0) {
				r1[r1<0] = -r1[r1<0]
			}
			if(rats[r[2]]<0) {
				r1[r1>0] = -r1[r1>0]
			}			
			so=(s[y]+lens[r[1]])-1
			d[s[y]:so] = r1
			rep=rbind(rep, cbind(chr,s[y],so, mean(r1), length(r1), var(r1)))
		}		
		data = NULL
		st=0
		for(x in 1:length(d)) {
			r = runif(1,1,length(spacing))
			so = spacing[r]
			st = st+so
			data=rbind(data,cbind(chr,st,st+60,d[x]))
		}	
		lis$data = data
		
		for(x in 1:length(rep[,1])) {
			rep[x,2] = data[rep[x,2],2]
			rep[x,3] = data[rep[x,3],3]
		}
		lis$rep = rep
	return(lis)
}
