`safe` <- function(set=NULL) {
	data = set$data
	if(is.null(data)) { data = read.table(set$files$file); }
	sa <- function(data, set=NULL) {
		n = unlist(strsplit(set$files$file, "/"))
		file = paste(set$files$odir, "/T", substr(n[length(n)], 0, nchar(n[length(n)])), sep="")
		write.table(data, file=file, sep="\t", row.names=F, col.names=F, quote=F)
		cur = getwd(); Jpath =  system.file("java", package="CNsolidate"); setwd(Jpath);
		command = paste("java -Xmx1600m safeDetect -f ", file, " -mr ",	set$safe$mR, " -np ", set$safe$mL, " -i ", set$files$idir, " -o ", set$files$odir, sep="")
		system(command); setwd(cur);
		na = paste(set$files$odir, "/T", substr(n[length(n)], 0, nchar(n[length(n)])-4), "_FinalReport_safe.txt", sep="");
		r = read.table(na);
		command = paste("rm ", na, " ", file, sep=""); system(command);
		return(r)
	}
	n = unlist(strsplit(set$files$file, "/"))
	name = paste(set$files$odir, "/", substr(n[length(n)], 0, nchar(n[length(n)])-4), "_FinalReport_safe.txt", sep="")
	r = rbind(sa(data[data[,4]>0,], set), sa(data[data[,4]<0,], set));
	write.table(r, file=name, sep="\t", row.names=F, col.names=F, quote=F)
	invisible(r)
}


#`safe.dec` <- function(file, idir, odir, mp=5, mr=0.2) {		
		
		#ifile = paste(idir, file, sep="")
#		cur = getwd()
#		Jpath =  system.file("java", package="CNsolidate")
#		setwd(Jpath)
#		command = paste("java -Xmx1600m safeDetect -f ", file, " -mr ", mr, " -np ", mp," -i ", idir, " -o ", odir, sep="")
#		system(command)
#		setwd(cur)
#	}

