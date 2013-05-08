`mapBreakc` <- function(set=NULL) {
	data = set$data; nrep=NULL;
	if(is.null(data)) { data = read.table(set$files$file); }	
	n = unlist(strsplit(set$files$file, "/"))
	repName = paste(set$files$odir, "/", substr(n[length(n)], 0, nchar(n[length(n)])-4), "_FinalReport.txt", sep="")
	if(length(count.fields(repName))>0) {
	report=read.table(repName);	
	u = unique(report[,1]); 
	for(x in 1:length(u)) {
		d=data[data[,1]==u[x],]; rep=report[report[,1]==u[x],]; pin = 1;
		v1=vector(); v2=vector(); rat = vector();
		
		for(y in 1:length(rep[,1])) {
			re=rep[y,]; ind1=re[6]; ind2=re[7];

			r <- .C("_MAP"
					,"dat" = as.double(d[,4])
					,"i1" = as.integer(ind1-1)
					,"i2" = as.integer(ind2-1)
					,"rat" = as.double(set$combine.settings$absRatio)
					,"size" = as.integer(length(d[,1]))
			)
			stoo = r$i2+1;
			if(stoo > length(d[,1])) { stoo=length(d[,1]); }
			v1[pin]=r$i1+1; v2[pin]=stoo; rat[pin]=r$rat; pin=pin+1;
		}
	v1[v1==0]=1; v2[v2>length(d[,1])]==length(d[,1]);	
	res=cbind(d[v1,1], d[v1,2], d[v2,3], rat, (v2-v1)+1, v1, v2);
		for(z in 1:length(res[,1])) {
			res[z,4] = mean(d[res[z,6]:res[z,7],4])
		}
		nrep=rbind(nrep,res);
	}
	if(length(nrep)>0) {
		if(length(nrep)==7) {
			nrep[1]=as.integer(nrep[1]); nrep[2]=as.integer(nrep[2]); nrep[3]=as.integer(nrep[3]);
			nrep=nrep[nrep[5]>set$combine.settings$minProbe-1 & abs(nrep[4])>set$combine.settings$absRatio]; 
			write.table(t(nrep), file=repName, sep="\t", row.names=F, col.names=F, quote=F)
		} else {
			nrep[,1]=as.integer(nrep[,1]); nrep[,2]=as.integer(nrep[,2]); nrep[,3]=as.integer(nrep[,3]); nrep=nrep[nrep[,5]>set$combine.settings$minProbe-1 & abs(nrep[,4])>set$combine.settings$absRatio,]; 
			if(length(nrep)==7) {
				write.table(t(nrep), file=repName, sep="\t", row.names=F, col.names=F, quote=F)
			} else {
				write.table(data.frame(as.integer(nrep[,1]), as.integer(nrep[,2]), as.integer(nrep[,3]), nrep[,-(1:3)]), file=repName, sep="\t", row.names=F, col.names=F, quote=F)
			}
		}
		
	} else {
		write.table(data.frame(as.integer(nrep[,1]), as.integer(nrep[,2]), as.integer(nrep[,3]), nrep[,-(1:3)]), file=repName, sep="\t", row.names=F, col.names=F, quote=F) 
		#write.table(nrep, file=repName, sep="\t", row.names=F, col.names=F, quote=F)
	}
	}
	return(nrep)
}


#`mapBreakc` <- function(file, odir, mr=0.3, mp=3) {	
#	n = unlist(strsplit(file, "/"))
#	repName = paste(odir, substr(n[length(n)], 0, nchar(n[length(n)])-4), "_FinalReport.txt", sep="")
#	data=read.table(file); report=read.table(repName);	
#	u = unique(report[,1]); nrep=NULL;
#	for(x in 1:length(u)) {
#		d=data[data[,1]==u[x],]; rep=report[report[,1]==u[x],]; pin = 1;
#		v1=vector(); v2=vector(); rat = vector();
#		
#		for(y in 1:length(rep[,1])) {
#			re=rep[y,]; ind1=re[6]; ind2=re[7];
#
#			r <- .C("_MAP"
#					,"dat" = as.double(d[,4])
#					,"i1" = as.integer(ind1)
#					,"i2" = as.integer(ind2)
#					,"rat" = as.double(mr)
#					,"size" = as.integer(length(d[,1]))
#			)
#			stoo = r$i2+1;
#			if(stoo > length(d[,1])) { stoo=length(d[,1]); }
#			v1[pin]=r$i1+1; v2[pin]=stoo; rat[pin]=r$rat; pin=pin+1;
#		}
#		
#	res=cbind(d[v1,1], d[v1,2], d[v2,3], rat, (v2-v1)+1, v1, v2);
#		for(z in 1:length(res[,1])) {
#			res[z,4] = mean(d[res[z,6]:res[z,7],4])
#		}
#		nrep=rbind(nrep,res);
#	}	
#	nrep=nrep[nrep[,5]>mp-1 & abs(nrep[,4])>mr,]; nrep[,1]=as.integer(nrep[,1]); nrep[,2]=as.integer(nrep[,2]); nrep[,3]=as.integer(nrep[,3]);
#	write.table(data.frame(as.integer(nrep[,1]),as.integer(nrep[,2]),as.integer(nrep[,3]), nrep[,-(1:3)]), file=repName, sep="\t", row.names=F, col.names=F, quote=F)	
#	return(nrep)
#}