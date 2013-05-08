`QC`<-function(set=NULL) {
	file=set$files$file; odir=set$files$odir;
	n = unlist(strsplit(file, "/"));
	name = paste(odir, "/", substr(n[length(n)], 0, nchar(n[length(n)])-4), "_QC.txt", sep="");
	rep = paste(odir, "/", substr(n[length(n)], 0, nchar(n[length(n)])-4), "_FinalReport.txt", sep="");
	report = read.table(rep, sep="\t"); d = read.table(file, sep="\t");
	r = quantile(abs(d[,4]), probs=0.68); dl = dLRs(d[,4]); w = ws(d[,4]);
	tot = length(report[report[,1]<23,1]); 
	vals = c(dl, r, w, tot)
	write.table(vals, file=name, sep="\t", row.names=F, col.names=F, quote=F)
}
