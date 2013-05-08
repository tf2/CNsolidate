`rId2` <-
function(d, chr, reso) {
	t = 0;
	chr1 = d[d[,1]==chr,]
	len = length(chr1[,1])
	fac = ceiling(len / reso)
	st = 1;
	sto = reso;
	for (a in 1:fac) {
		dat = chr1[st:sto,4]
		dat = as.numeric(dat)
		if (length(dat)<100) {
			zer = vector(length=100)
			for (z in 1:length(zer)) {
					zer[z] = 0
			}
				 temp = c(dat,zer);
				 te = bcp(temp)$posterior.prob
				 pin = length(te)-100
			t = c(t, te[1:pin])
		}
		else if (length(chr1[,1]) < reso) {
        	dat = chr1[1:len, 4]
        	t = c(t, smugFun(dat))
        	
		} else {
			t =  c(t,bcp(dat)$posterior.prob);
		}
		st = st+reso;
		sto = sto+reso;
		
		if (sto>len) {
			sto = len;
		}
	}
	t=t[-(1)]
	return(cbind(chr1[,1:4],t))
}

