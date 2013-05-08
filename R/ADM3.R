`ADM3` <- function(set=NULL) {
	data = set$data
	if(is.null(data)) { data = read.table(set$files$file); }
	report=NULL; u=unique(data[,1]);
	for(x in 1:length(u)) {
		print(u[x])
		d=data[data[,1]==u[x],]; d=ADM_score(d); dd=ADM_segment(d,set$ADM3$t,set$ADM3$nd); ddd=ADM_feature(dd, set$ADM3$np);
		report=rbind(report,ddd);
	}
	report=report[abs(report[,4])>set$ADM3$mR & report[,7]>set$ADM3$np,];
	l=list(); l$data=data; l$report=report;
	n = unlist(strsplit(set$files$file, "/"))
	name = paste(set$files$odir, "/", substr(n[length(n)], 0, nchar(n[length(n)])-4), "_FinalReport_ADM3.txt", sep="")
	n=l$report; if(length(n)==8) { n = rbind(n,n); }
	write.table(n, file=name, sep="\t", row.names=F, col.names=F, quote=F)
	invisible(l)
}


#`ADM3` <- function(data, file, odir, t=5, nd=2, mr=0.4, np=2) {	
#	report=NULL; u=unique(data[,1]);
#	for(x in 1:length(u)) {
#		print(x)
#		d=data[data[,1]==u[x],]; d=ADM_score(d); dd=ADM_segment(d,t,nd); ddd=ADM_feature(dd, np);
#		report=rbind(report,ddd);
#	}
#	report=report[abs(report[,4])>mr & report[,7]>np,];
#	l=list(); l$data=data; l$report=report;
#	n = unlist(strsplit(file, "/"))
#	name = paste(odir, "P", substr(n[length(n)], 0, nchar(n[length(n)])-4), "_FinalReport_ADM3.txt", sep="")
#	p=l$report[l$report[,4]>0,]; if(length(p)==8) { p = rbind(p,p); }
#	write.table(p, file=name, sep="\t", row.names=F, col.names=F, quote=F)
#	name = paste(odir, "N", substr(n[length(n)], 0, nchar(n[length(n)])-4), "_FinalReport_ADM3.txt", sep="")
#	n=l$report[l$report[,4]<=0,]; if(length(n)==8) { n = rbind(n,n); }
#	write.table(n, file=name, sep="\t", row.names=F, col.names=F, quote=F)
#	return(l)
#}
