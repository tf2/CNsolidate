`split_array` <- function(file, odir) {
	m=read.table(file, sep="\t"); d1=NULL; d2=NULL;
	m=m[order(m[,1], m[,2], m[,3]),]; u = unique(m[,1]);
	for(x in 1:length(u)) {
		d = m[m[,1]==u[x],];
		if(length(d[,1])<=5000) {
			d2= rbind(d2,d);
		} else {
			d1=rbind(d1,d);
		}
	}
	write.table(d1, file=file, sep="\t", row.names=F, col.names=F, quote=F)
	n = unlist(strsplit(file, "/"));
	name = paste(odir, "/", substr(n[length(n)], 0, nchar(n[length(n)])-4), "_duplicatedProbes.txt", sep="");
	write.table(d2, file=name, sep="\t", row.names=F, col.names=F, quote=F)
}