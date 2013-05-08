`chr_m`<-function(file, odir) {
	n = unlist(strsplit(file, "/"));
	name = paste(odir, "/", substr(n[length(n)], 0, nchar(n[length(n)])-4), "_Chr_Medians.txt", sep="");
	d = read.table(file, sep="\t"); d[,4]=d[,4]-median(d[d[,1]<23,4], na.rm=T);
	vals = vector(); u = unique(d[,1])
	for(x in 1:length(u)) {
		vals[x] = median(d[d[,1]==u[x],4], na.rm=T)
	}
	write.table(cbind(u, vals), file=name, sep="\t", row.names=F, col.names=F, quote=F)
	return(vals)
}
