`fsegment` <- function(set=NULL) {
	data = set$data
	if(is.null(data)) { data = read.table(set$files$file); }
		fseg <- function(data, set) {
			uni = unique(data[,1]); rep = NULL;	
			for(x in 1:length(uni)) {
				dat = data[data[,1]==uni[x],]; dat = dat[order(dat[,2],dat[,3]),];
				v = lowcost.matrix.rep(dat, set$fsegment$sL, set$fsegment$cpP); rep = rbind(rep, v)
			}
			if(!is.null(rep)) {
				#rep[,1] = as.integer(rep[,1]); rep[,2] = as.integer(rep[,2]); rep[,3] = as.integer(rep[,3])
				rep = data.frame(as.integer(rep[,1]), as.integer(rep[,2]), as.integer(rep[,3]), rep[,4], rep[,5])
			}
		invisible(rep)
		}	
		n = unlist(strsplit(set$files$file, "/")); name = paste(set$files$odir, "/", substr(n[length(n)], 0, nchar(n[length(n)])-4), sep = "")
		rep = NULL;
		if(set$algorithms.sep$fsegment==1) {
			r1 = fseg(data[data[,4]>0,], set); r2 = fseg(data[data[,4]<0,], set); rep = rbind(r1,r2);
		} else {
			rep = fseg(data, set);
		}
		if(!is.null(rep)) {
			rep=rep[abs(rep[,4])>set$fsegment$mR & rep[,5]>=set$fsegment$mL,]
			rep=rep[order(rep[,1], rep[,2], rep[,3]),];
		}
		write.table(rep, file=paste(name, "_FinalReport_fsegment.txt", sep=""), sep="\t", quote=F, col.names=F, row.names=F)
invisible(rep)
}



#`fsegment` <- function(data, file, idir, odir, segLen, cpPex, minL=3, mr=0.3) {
#		
#		uni = unique(data[,1])
#		rep = NULL
#		
#		for(x in 1:length(uni)) {
#			my = data[data[,1]==uni[x],]
#			my = my[order(my[,2],my[,3]),]
#			v = lowcost.matrix.rep(my, segLen, cpPex)
#			rep = rbind(rep, v)
#		}
#
#		n = unlist(strsplit(file, "/"))
#		name = paste(odir, substr(n[length(n)], 0, nchar(n[length(n)])-4), sep = "")
#		rep=rep[abs(rep[,4])>mr & rep[,5]>=minL,]
#		rep[,1] = as.integer(rep[,1])
#		rep[,2] = as.integer(rep[,2])
#		rep[,3] = as.integer(rep[,3])
#		write.table(rep, file=paste(name, "_FinalReport_fsegment.txt", sep=""), sep="\t", quote=F, col.names=F, row.names=F)
#	
#return(rep)
#}
