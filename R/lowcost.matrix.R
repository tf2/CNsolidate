 `lowcost.matrix` <- function(dat, segLen, cpPex) {
 
 	s = seq(1, length(dat), by=segLen)
 	so = segLen
 	cpPex = (segLen/100)*cpPex
 	vec = vector(length=length(dat))

 	for(x in 1:length(s)) {
 		st = s[x]
 		tst = st
 		ss = 1
 		if (so > length(dat)) {
 			so = length(dat)
 			tlen = so-st
 			cpPex = (tlen/100)*cpPex
 			if(cpPex < 1) {
 				cpPex=1
 			}
 		}
 		m = dat[st:so]
 		res = findsegments(m, maxcp=cpPex, maxk=length(m), verbose=0)
 		p = res$th[length(res$th[,1]),]

 			for(y in 1:length(p)) {
 				soo = p[y] -1				
 				tso = (st+soo)-1
 				mu = mean(m[ss:soo])
 				ss = p[y] 				 				
 				vec[tst:tso] = mu
 				tst = st+p[y]
 			}
 		so=so+segLen
 	}
 	print(so)
 	return(vec)
 }
 
 
 
 