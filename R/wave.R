`wave` <- function (set=NULL) {
	m = read.table(set$files$file);
	cur = getwd(); Jpath = system.file("java", package = "CNsolidate"); setwd(Jpath);
    wfile1 = paste(set$files$odir, "/waveTemp.txt", sep=""); wfile2 = paste(set$files$odir, "/waveTemp1.txt", sep="")
    m = m[order(m[, 1], m[, 2], m[, 3]), ]; u = unique(m[, 1]); r = NULL;
    log2_dye_self=m[,4]
    	for (x in 1:length(u)) {
       		d = m[m[, 1] == u[x], ]; me = median(d[, 4]); d[, 4] = d[, 4] - me;
       		if (length(d[, 1]) > 5000) {
       			l = round(length(d[,4])/10)
       			ll = round(length(d[,4])/500)
       			if(!l%% 2) { l=l+1; }
       			if(!ll%% 2) { ll=ll+1; }
       			rem = runmed(d[,4], ll)-runmed(d[,4], l)
       			rem=rem[abs(rem)<quantile(abs(rem), probs=0.68)]
       			wf = (max(rem)-min(rem))*14
       			print(wf)
       			rr = rem>0; n = 0; dist = 0; check=1;
				for (x in 1:length(rr)) {
					if(rr[x]) {
						check=1; dist=dist+1;
					} else {
						if(check) {
							n=n+1; check=0;
						}
					}
				}
           		write.table(d, file = wfile1, sep = "\t", row.names = F, col.names = F, quote = F)
           		command = paste("java -Xmx1600m Wave -f ", wfile1, " -fa ", wf, " -si ", (dist/n)*2, " > ", wfile2, sep = "")
           		system(command); w = read.table(wfile2); d[, 4] = w[, 3] + me;
           		r = rbind(r, d)
        	} else {
            	r = rbind(r, d)
        	}
    	}
    command = paste("rm ", wfile1, " ", wfile2, sep=""); system(command); setwd(cur);
	r=cbind(r[,1:5], log2_dye_self, r[,-(1:5)])
	write.table(r, file=set$files$file, sep="\t", row.names=F, col.names=F, quote=F)
invisible(r)
}

`wavef` <- function(set=NULL) {
	m = read.table(set$files$file);
	cur = getwd(); Jpath = system.file("java", package = "CNsolidate"); setwd(Jpath);
	wfile1 = paste(set$files$odir, "/waveTemp.txt", sep=""); wfile2 = paste(set$files$odir, "/waveTemp1.txt", sep="")
	m=m[order(m[,1], m[,2], m[,3]),]; u = unique(m[,1]); r = NULL
	log2_dye_self=m[,4]
	for(x in 1:length(u)) {
		d = m[m[,1]==u[x],]; me = median(d[,4]); d[,4]=d[,4]-me
		if(length(d[,1])>5000) {
			write.table(d, file=wfile1, sep="\t", row.names=F, col.names=F, quote=F)
			command =paste("java -Xmx1600m FixedWave -f ",wfile1, " -fa ", set$norm.settings$wFac, " > ", wfile2, sep="" ); system(command)
			w=read.table(wfile2); d[,4]=w[,3]+me; r=rbind(r,d);
		} else {
			r=rbind(r,d);
		}
	}
	command = paste("rm ", wfile1, " ", wfile2, sep=""); system(command); setwd(cur);
	r=cbind(r[,1:5], log2_dye_self, r[,-(1:5)])
	write.table(r, file=set$files$file, sep="\t", row.names=F, col.names=F, quote=F)
}

#`wave` <- function (set=NULL) {
#	m = read.table(set$files$file);
#	cur = getwd(); Jpath = system.file("java", package = "CNsolidate"); setwd(Jpath);
#   wfile1 = paste(set$files$odir, "/waveTemp.txt", sep=""); wfile2 = paste(set$files$odir, "/waveTemp1.txt", sep="")
#    m = m[order(m[, 1], m[, 2], m[, 3]), ]; u = unique(m[, 1]); r = NULL;
#    log2_dye_self=m[,4]
#    	for (x in 1:length(u)) {
#       		d = m[m[, 1] == u[x], ]; me = median(d[, 4]); d[, 4] = d[, 4] - me;
#       		if (length(d[, 1]) > 5000) {
#       			l = round(length(m[m[, 1] == u[x],4])/10)
#       			ll = round(length(m[m[, 1] == u[x],4])/500)
#       			if(l%% 2) { l=l+1; }
#       			if(ll%% 2) { ll=ll+1; }
#       			rem = runmed(m[m[, 1] == u[x],4], ll)-runmed(m[m[, 1] == u[x],4], l)
#       			rem=rem[abs(rem)<quantile(abs(rem), probs=0.68)]
#       			wf = (max(rem)-min(rem))*14
#       			print(wf)
#           		write.table(d, file = wfile1, sep = "\t", row.names = F, col.names = F, quote = F)
#           		command = paste("java -Xmx1600m Wave -f ", wfile1, " -fa ", wf, " > ", wfile2, sep = "")
#           		system(command); w = read.table(wfile2); d[, 4] = w[, 3] + me;
#           		r = rbind(r, d)
#        	} else {
#            	r = rbind(r, d)
#        	}
#    	}
#    command = paste("rm ", wfile1, " ", wfile2, sep=""); system(command); setwd(cur);
#	r=cbind(r[,1:5], log2_dye_self, r[,-(1:5)])
#	write.table(r, file=set$files$file, sep="\t", row.names=F, col.names=F, quote=F)
#invisible(r)
#}

#`wave` <- function(set=NULL) {
#	m = read.table(set$files$file);
#	cur = getwd(); Jpath = system.file("java", package = "CNsolidate"); setwd(Jpath);
#	wfile1 = paste(set$files$odir, "/waveTemp.txt", sep=""); wfile2 = paste(set$files$odir, "/waveTemp1.txt", sep="")
#	m=m[order(m[,1], m[,2], m[,3]),]; u = unique(m[,1]); r = NULL
#	log2_dye_self=m[,4]
#	for(x in 1:length(u)) {
#		d = m[m[,1]==u[x],]; me = median(d[,4]); d[,4]=d[,4]-me
#		if(length(d[,1])>5000) {
#			write.table(d, file=wfile1, sep="\t", row.names=F, col.names=F, quote=F)
#			command =paste("java -Xmx1600m Wave -f ",wfile1, " -fa ", set$norm.settings$wFac, " > ", wfile2, sep="" ); system(command)
#			w=read.table(wfile2); d[,4]=w[,3]+me; r=rbind(r,d);
#		} else {
#			r=rbind(r,d);
#		}
#	}
#	command = paste("rm ", wfile1, " ", wfile2, sep=""); system(command); setwd(cur);
#	r=cbind(r[,1:5], log2_dye_self, r[,-(1:5)])
#	write.table(r, file=set$files$file, sep="\t", row.names=F, col.names=F, quote=F)
#}


#`wave` <- function(file, idir, odir, fa=0.5) {
#	m = read.table(file);
#	cur = getwd()
#   Jpath = system.file("java", package = "CNsolidate")
#   setwd(Jpath)
#	wfile1 = paste(odir, "/waveTemp.txt", sep=""); wfile2 = paste(odir, "/waveTemp1.txt", sep="")
#	m=m[order(m[,1], m[,2], m[,3]),]; u = unique(m[,1]); r = NULL
#	for(x in 1:length(u)) {
#		d = m[m[,1]==u[x],]; me = median(d[,4]); d[,4]=d[,4]-me
#		if(length(d[,1])>5000) {
#			write.table(d, file=wfile1, sep="\t", row.names=F, col.names=F, quote=F)
#			command =paste("java -Xmx1600m Wave -f ",wfile1, " -fa ", fa, " > ", wfile2, sep="" ); system(command)
#			w=read.table(wfile2);
#			d[,4]=w[,3]+me; r=rbind(r,d);
#		} else {
#			r=rbind(r,d);
#		}
#		
#	}
#	command = paste("rm ", wfile1, " ", wfile2, sep=""); system(command);
#	setwd(cur);
#	write.table(r, file=file, sep="\t", row.names=F, col.names=F, quote=F)
#}
