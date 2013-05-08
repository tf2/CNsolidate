`gap`<-function(set=NULL) {
	n = unlist(strsplit(set$files$file, "/"))
	rep = paste(set$files$odir, "/", substr(n[length(n)], 0, nchar(n[length(n)])-4), "_FinalReport.txt", sep="")
	if(length(count.fields(rep))>0) {
		cur = getwd(); Jpath =  system.file("java", package="CNsolidate"); setwd(Jpath);
		command = paste("java -Xmx1600m Gap -f ", set$files$file, " -r ", rep, " -g ", set$files$gapfile, sep="")
		system(command); setwd(cur); r=read.table(rep); r=r[r[,5]>set$combine.settings$minProbe-1 & abs(r[,4])>set$combine.settings$absRatio,];
		write.table(data.frame(as.integer(r[,1]), as.integer(r[,2]), as.integer(r[,3]), r[,-(1:3)]), file=rep, sep="\t", row.names=F, col.names=F, quote=F)
	}
}


#`gap`<-function(file, idir, odir, mr=0.3, mp=3, gold = "../extdata/centromeres_hg19.txt") {
#		n = unlist(strsplit(file, "/"))
#		rep = paste(odir, substr(n[length(n)], 0, nchar(n[length(n)])-4), "_FinalReport.txt", sep="")
#		cur = getwd()
#		Jpath =  system.file("java", package="CNsolidate")
#		setwd(Jpath)
#		command = paste("java -Xmx1600m Gap -f ", file, " -r ", rep, " -g ", gold, sep="")
#		system(command)
#		setwd(cur)
#		r=read.table(rep); r=r[r[,5]>mp-1 & abs(r[,4])>mr,];
#		write.table(data.frame(as.integer(r[,1]), as.integer(r[,2]), as.integer(r[,3]), r[,-(1:3)]), file=rep, sep="\t", row.names=F, col.names=F, quote=F)
#}