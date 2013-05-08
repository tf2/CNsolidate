`SW.array` <- function(data, reso=1000) {
		
	u = unique(data[,1]) 
	res = NULL
	for(i in 1:length(u)) {
	
		print(u[i])
		chr1 = data[data[,1]==u[i],]
		len = length(chr1[,1])
		fac = ceiling(len / reso)
		st = 1;
		sto = reso;
		
		for (a in 1:fac) {

			dat = chr1[st:sto,]
			len1 = length(dat[,1])
			dat[,1] = 1
			
			if (len1<reso) {
			
				zer = matrix(nrow=reso, ncol=length(dat[1,]))
				
				for (z in 1:length(zer[,1])) {
					for(zz in 1:length(zer[1,])) {
						zer[z,zz] = mean(dat[,zz])
					}
				}
			
				dat1 = rbind(dat, zer)
				te = hmm.run.func(dat1)
				te = te[[1]]
				te = te[[1]]
				te = cbind(te[1:len1,3],te[1:len1,6])		
				st = st+reso;
				sto = sto+reso;
			
			} else {
			
			te = hmm.run.func(dat)
			te = te[[1]]
			te = te[[1]]
			te = cbind(te[,3],te[,6])
			st = st+reso;
			sto = sto+reso;
			
			}
			
			if (sto>len) {
				sto = len;
			}
			res = rbind(res, cbind(u[i], dat[,4:6], te)) 
		}

	}
	return(res)
}