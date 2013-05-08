`rcncp` <- function(set=NULL) {
		data = set$data
		if(is.null(data)) { data = read.table(set$files$file); }	
		rc <- function(data, set) {
			l = cncp(data,set$rcncp$fac)
			fea = cncnps(l, set$rcncp$a, set$rcncp$n, set$rcncp$prob, set$rcncp$dist)
			n = unlist(strsplit(set$files$file, "/"))
			name = paste(set$files$odir, substr(n[length(n)], 0, nchar(n[length(n)])-4), "_FinalReport_rcncp.txt", sep="")
			if(!is.null(dim(fea))) {
				fea=fea[abs(fea[,4])>set$rcncp$mR & abs(fea[,5])>set$rcncp$a & abs(fea[,6])>set$rcncp$prob,]
			}
			if(!is.null(dim(fea))) {
				fea= data.frame(as.integer(fea[,1]), as.integer(fea[,2]), as.integer(fea[,3]), fea[,-(1:3)])
				write.table(fea, file=name, sep="\t", quote=F, col.names=F, row.names=F)
			} else {
				fea =""
				write.table(fea, file=name, sep="\t", quote=F, col.names=F, row.names=F)
			}
		}
		n = unlist(strsplit(set$files$file, "/"))
		name = paste(set$files$odir, substr(n[length(n)], 0, nchar(n[length(n)])-4), "_FinalReport_rcncp.txt", sep="")
		r=NULL;
		if(set$algorithms.sep$rcncp==1) {
			r1 = NULL; r2 = NULL;
			rc(data[data[,4]>0,], set); 
			if(length(count.fields(name))>0) {
				r1 = read.table(name);
			}
			rc(data[data[,4]<0,], set); 
			if(length(count.fields(name))>0) {
				r2 = read.table(name);
			}
			r = rbind(r1,r2); 
		} else {
			rc(data, set); 
			if(length(count.fields(name))>0) {
				r = read.table(name);
			}
		}
		if(length(r[,1])>0) {
			r=r[order(r[,1],r[,2],r[,3]),];
		}
		write.table(r, file=name, sep="\t", quote=F, col.names=F, row.names=F)
}



#`rcncp` <- function(data, file, idir, odir, fac = 3, a = 3, n = 2, absc = 0.3, prob=0.4, dist=5000) {
#	
#		l = cncp(data,fac)
#		fea = cncnps(l, a, n, prob, dist)
#		if(!is.null(dim(fea))) {
#			fea=fea[abs(fea[,4])>absc & abs(fea[,5])>a & abs(fea[,6])>prob,]
#		}
#		if(!is.null(dim(fea))) {
#			fea[,1] = as.integer(fea[,1])
#			fea[,2] = as.integer(fea[,2])
#			fea[,3] = as.integer(fea[,3])
#			n = unlist(strsplit(file, "/"))
#			name = paste(odir, substr(n[length(n)], 0, nchar(n[length(n)])-4), "_FinalReport_CNCP.txt", sep="")
#			write.table(fea, file=name, sep="\t", quote=F, col.names=F, row.names=F)
#		} else {
#			fea =""
#			n = unlist(strsplit(file, "/"))
#			name = paste(odir, substr(n[length(n)], 0, nchar(n[length(n)])-4), "_FinalReport_CNCP.txt", sep="")
#			write.table(fea, file=name, sep="\t", quote=F, col.names=F, row.names=F)
#			
#		}
#}