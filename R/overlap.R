`overlap` <- function(set=NULL) {		
		n = unlist(strsplit(set$files$file, "/"))
		file1 = paste(set$files$odir, "/", substr(n[length(n)], 0, nchar(n[length(n)])-4), "_Allfeatures.txt", sep="")
		file2 = paste(set$files$odir, "/", substr(n[length(n)], 0, nchar(n[length(n)])-4), "_FinalReport.txt", sep="")
		if(length(count.fields(file1))>0 & length(count.fields(file2))>0) {
		cur = getwd(); Jpath =  system.file("java", package="CNsolidate"); setwd(Jpath);
		command = paste("java -Xmx1600m Overlap -f1 ", file1, " -f2 ", file2, " -g ", set$files$gold, sep="")
		system(command); setwd(cur); r=read.table(file2); 
		write.table(data.frame(as.integer(r[,1]), as.integer(r[,2]), as.integer(r[,3]), r[,-(1:3)]), file=file2, sep="\t", row.names=F, col.names=F, quote=F)
		}
}

#`overlap` <- function(file, idir, odir, gold = "../extdata/all_sets_merged_hg19.txt") {		
#		n = unlist(strsplit(file, "/"))
#		file1 = paste(odir, substr(n[length(n)], 0, nchar(n[length(n)])-4), "_Allfeatures.txt", sep="")
#		file2 = paste(odir, substr(n[length(n)], 0, nchar(n[length(n)])-4), "_FinalReport.txt", sep="")
#		cur = getwd()
#		Jpath =  system.file("java", package="CNsolidate")
#		setwd(Jpath)
#		command = paste("java -Xmx1600m Overlap -f1 ", file1, " -f2 ", file2, " -g ", gold, sep="")
#		system(command)
#		setwd(cur)
#		r=read.table(file2); 
#		write.table(data.frame(as.integer(r[,1]), as.integer(r[,2]), as.integer(r[,3]), r[,-(1:3)]), file=file2, sep="\t", row.names=F, col.names=F, quote=F)
#}
