 `lowcost.matrix.rep` <- function(dat, segLen, cpPex) {
 
 	s = seq(1, length(dat[,1]), by=segLen)
 	so = segLen
 	cpPex = (segLen/100)*cpPex
 	starts = NULL
 	stops = NULL
 	mens = NULL
 	lens = NULL

 	for(x in 1:length(s)) {
 		st = s[x]
 		tst = st
 		ss = 1
 		if (so > length(dat[,1])) {
 			so = length(dat[,1])
 			tlen = so-st
 			cpPex = (tlen/100)*cpPex
 			if(cpPex < 1) {
 				cpPex=1
 			}
 		}

 		m = dat[st:so,4]
 		if(length(m) > 1) { 
 		res = findsegments(m, maxcp=cpPex, maxk=length(m), verbose=0)
 		p = res$th[length(res$th[,1]),]

 			for(y in 1:length(p)) {
 				soo = p[y] -1				
 				tso = (st+soo)-1
 				mu = mean(m[ss:soo])
 				lens = c(lens, length(m[ss:soo]))
 				ss = p[y]
 				if (tst> length(dat[,1])) {
 					tst=length(dat[,1])-1
 				}
 				if (tso> length(dat[,1])) {
 					tso=length(dat[,1])
 				}
 				starts = c(starts, dat[tst,2])
 				stops = c(stops, dat[tso,3])
 				mens = c(mens, mu)
 				tst = st+p[y]
 			}
 			so=so+segLen
 		} 
 		
 	}
 	print(so)
 	return(cbind(dat[1,1], starts, stops, mens, lens))
 }
 
 
 
 