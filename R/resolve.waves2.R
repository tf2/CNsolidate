`resolve.waves2` <- function (set=NULL) {
	rr = get.wave.form(set)
	m = read.table(set$files$file);
	cur = getwd(); Jpath = system.file("java", package = "CNsolidate"); setwd(Jpath);
	wfile1 = paste(set$files$odir, "/waveTemp.txt", sep=""); wfile2 = paste(set$files$odir, "/waveTemp1.txt", sep="")
    m = m[order(m[, 1], m[, 2], m[, 3]), ]; u = unique(m[, 1]); r = NULL;
    log2_dye_self=m[,4]	
    	for (x in 1:length(u)) {
       		d = m[m[, 1] == u[x], ]; me = median(d[, 4]); d[, 4] = d[, 4] - me;
       		if (length(d[, 1]) > 5000) {
       			l = round(length(m[m[, 1] == u[x],4])/10)
       			if(l%% 2) { l=l+1; }
       			ww = m[m[, 1] == u[x],4] - runmed(m[m[, 1] == u[x],4], l)
       			wf= quantile(abs(ww), 0.68)*(quantile(abs(ww), 0.68)*45)
       			print(wf)
           		write.table(d, file = wfile1, sep = "\t", row.names = F, col.names = F, quote = F)
           		command = paste("java -Xmx1600m Wave -f ", wfile1, " -fa ", wf, " > ", wfile2, sep = "")
           		system(command); w = read.table(wfile2); d[, 4] = w[, 3] + me;
           		r = rbind(r, d)
        	} else {
            	r = rbind(r, d)
        	}
    	}
    command = paste("rm ", wfile1, " ", wfile2, sep=""); system(command); setwd(cur);
	r=cbind(r[,1:5], log2_dye_self, r[,-(1:5)])
	write.table(r, file=set$files$file, sep="\t", row.names=F, col.names=F, quote=F)
}
