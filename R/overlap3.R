`overlap3` <- function(set=NULL) {		
		n = unlist(strsplit(set$files$file, "/"))
		file = paste(set$files$odir, "/", substr(n[length(n)], 0, nchar(n[length(n)])-4), "_FinalReport.txt", sep="")
		cur = getwd(); Jpath =  system.file("java", package="CNsolidate"); setwd(Jpath);
		gold= paste(substr(set$files$gold, 1, nchar(set$files$gold)-4), "_deletion.txt", sep="")
		command = paste("java -Xmx1600m Overlap3 -f ", file, " -o ", file, " -g ", gold, sep="")
		system(command); 
		gold= paste(substr(set$files$gold, 1, nchar(set$files$gold)-4), "_duplication.txt", sep="")
		command = paste("java -Xmx1600m Overlap3 -f ", file, " -o ", file, " -g ", gold, sep="")
		system(command);
		gold= paste(substr(set$files$gold, 1, nchar(set$files$gold)-4), "_complex.txt", sep="")
		command = paste("java -Xmx1600m Overlap3 -f ", file, " -o ", file, " -g ", gold, sep="")
		system(command); 
		setwd(cur); r=read.table(file); 
		write.table(data.frame(as.integer(r[,1]), as.integer(r[,2]), as.integer(r[,3]), r[,-(1:3)]), file=file, sep="\t", row.names=F, col.names=F, quote=F)
}