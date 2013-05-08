`walk` <- function(set=NULL) {
    cur = getwd(); Jpath = system.file("java", package = "CNsolidate"); setwd(Jpath)
    data = set$data
	if(is.null(data)) { data = read.table(set$files$file); }
    n = unlist(strsplit(set$files$file, "/"))
   	f =  paste(set$files$odir, "/", substr(n[length(n)], 0, nchar(n[length(n)]) - 4), "_variancewalking.temp", sep = "")
   	name = paste(set$files$odir, "/", substr(n[length(n)], 0, nchar(n[length(n)]) - 4), "_variancewalking_walk.temp", sep = "")
    out = paste(set$files$odir, "/",substr(n[length(n)], 0, nchar(n[length(n)]) - 4), "_FinalReport_walk.txt", sep = "");
    u = unique(data[,1]); dd=NULL;
    for(x in 1:length(u)) {
    	if(u[x] < 23) {	
    		write.table(data[data[,1]==u[x],], file=f, row.names=F, col.names=F, quote=F, sep="\t")
    		command = paste("java -Xmx1600m ControlT -f ", f, " -mr ", set$walk$mR, " -i ", set$files$idir, " -o ", set$files$odir, sep = ""); system(command)
    		command = paste("java -Xmx1600m Clean -f ", name, " -mp ", set$walk$mL, " -mr ", set$walk$mR, " -pc ", set$walk$pc, " -pt ", set$walk$pT, " -gold ", set$files$gold, " -ty walk", sep = ""); system(command) 
    		if(length(count.fields(name))>0) {
    			dd=rbind(dd, read.table(name))
    		}
    	}
    }
    setwd(cur); command = paste("rm", f, name, sep=" "); system(command)
    write.table(dd, file=out, row.names=F, col.names=F, sep="\t", quote=F)
}