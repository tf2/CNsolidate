`get.wave.form` <- function (set=NULL) {
	m = read.table(set$files$file);
	cur = getwd(); Jpath = system.file("java", package = "CNsolidate"); setwd(Jpath);
	wfile1 = paste(set$files$odir, "/waveTemp.txt", sep=""); wfile2 = paste(set$files$odir, "/waveTemp1.txt", sep="")
    m = m[order(m[, 1], m[, 2], m[, 3]), ]; u = unique(m[, 1]); r = NULL;
    	for (x in 1:length(u)) {
       		d = m[m[, 1] == u[x], ]; me = median(d[, 4]); d[, 4] = d[, 4] - me;
       		if (length(d[, 1]) > 5000) {
       			wf=1.2
           		write.table(d, file = wfile1, sep = "\t", row.names = F, col.names = F, quote = F)
           		command = paste("java -Xmx1600m Wave -f ", wfile1, " -fa ", wf, " > ", wfile2, sep = "")
           		system(command); w = read.table(wfile2); d[, 4] =w[,2];
           		r = rbind(r, d)
        	} else {
            	r = rbind(r, d)
        	}
    	}
    command = paste("rm ", wfile1, " ", wfile2, sep = ""); system(command);
    setwd(cur);
return(r)
}
