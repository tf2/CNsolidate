`pn.signs` <-
function(data, file, idir, odir) {		

  		n = unlist(strsplit(file, "/"))	 
		pfin = paste(idir, "P", n[length(n)], sep="")
		nfin = paste(idir, "N", n[length(n)], sep="")
		
		write.table(data[data[,4]<=0,], file=nfin, sep="\t", col.names=F, row.names=F, quote=F)	
		write.table(data[data[,4]>0,], file=pfin, sep="\t", col.names=F, row.names=F, quote=F)	
	
	return(c(pfin, nfin))
}

