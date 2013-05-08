`SMA` <- function(set=NULL) {	
	data = set$data
	if(is.null(data)) { data = read.table(set$files$file); }
	rep=NULL; rr=NULL; u=unique(data[,1]);
	if(length(data[1,])<5) { data=cbind(data,1); }
	for(y in 1:length(u)) {
		print(u[y]); d=data[data[,1]==u[y],]; re=NULL;		  
		  res <- .C("_F",
			"d" =  as.double( r <- .C("_P"
							,"dat" = as.double(d[,4])
							,"prob" = as.double(d[,5])
		 					,"threshold" = as.double(set$SMA$s)
		 					,"size" = as.integer(length(d[,4]))
		 					,"PACKAGE" = "CNsolidate")$dat*(dLRs(d[,4])*set$SMA$fact)
		 					)
			,"st"= as.double(rep(0,length(d[,1])))
			, "so"= as.double(rep(0,length(d[,1])))
			,"threshold" = as.double(set$SMA$s)
			,"size" = as.integer(length(r))
			,"PACKAGE" = "CNsolidate"
			); ind1 = res$st[res$st!=0]; ind2 = res$so[res$so!=0];
			if(length(ind1)>0 & length(ind1) == length(ind2)) {
				for(x in 1:length(ind1)) {
					if(ind1[x]>0 & ind2[x]>0 & ind1[x]<=length(d[,1]) & ind2[x]<=length(d[,1]) & ind1[x]<=ind2[x]) {
						ind=mapBreak(d, ind1[x], ind2[x], set$SMA$mR);
						re = rbind(re, cbind(d[ind[1],1],  d[ind[1],2],  d[ind[2],3], mean(d[ind[1]:ind[2],4]), mean(r[ind[1]:ind[2]]), (ind[2]-ind[1])+1, ind[1], ind[2]))
					}
				};	
				rep=rbind(rep, re[abs(re[,4])>getCut(dLRs(d[,4]),set$SMA$mR) & re[,6]>set$SMA$np,]);
			}; rr=c(rr,res$d)
	}
	l=list(); l$data=cbind(data,rr); l$report=rep; set$data=l$data;
	n = unlist(strsplit(set$files$file, "/"))
	name = paste(set$files$odir, "/", substr(n[length(n)], 0, nchar(n[length(n)])-4), "_FinalReport_SMA.txt", sep="")
	n=l$report; if(length(n)==8) { n = rbind(n,n); }
	write.table(n, file=name, sep="\t", row.names=F, col.names=F, quote=F)
	invisible(set)
}

#`ADM2` <- function(data, file, odir, t=6, np=3, mr=0.3, fact=3) {		
#	rep=NULL; rr=NULL; u=unique(data[,1]);
#	if(length(data[1,])<5) { data=cbind(data,1); }
#	for(y in 1:length(u)) {
#		print(y); d=data[data[,1]==u[y],]; re=NULL;		  
#		  res <- .C("_F",
#			"d" =  as.double( r <- .C("_P"
#							,"dat" = as.double(d[,4])
#							,"prob" = as.double(d[,5])
#		 					,"threshold" = as.double(t)
#		 					,"size" = as.integer(length(d[,4]))
#		 					,"PACKAGE" = "CNsolidate")$dat*(dLRs(d[,4])*fact)
#		 					)
#			,"st"= as.double(rep(0,length(d[,1])))
#			, "so"= as.double(rep(0,length(d[,1])))
#			,"threshold" = as.double(t)
#			,"size" = as.integer(length(r))
#			,"PACKAGE" = "CNsolidate"
#			); ind1 = res$st[res$st!=0]; ind2 = res$so[res$so!=0];
#			if(length(ind1)>0 & length(ind1) == length(ind2)) {
#				for(x in 1:length(ind1)) {
#					if(ind1[x]>0 & ind2[x]>0 & ind1[x]<=length(d[,1]) & ind2[x]<=length(d[,1]) & ind1[x]<=ind2[x]) {
#						ind=mapBreak(d, ind1[x], ind2[x], mr);
#						re = rbind(re, cbind(d[ind[1],1],  d[ind[1],2],  d[ind[2],3], mean(d[ind[1]:ind[2],4]), mean(r[ind[1]:ind[2]]), (ind[2]-ind[1])+1, ind[1], ind[2]))
#					}
#				};	
#				rep=rbind(rep, re[abs(re[,4])>getCut(dLRs(d[,4]),mr) & re[,6]>np,]);
#			}; rr=c(rr,res$d)
#	}
#	l=list(); l$data=cbind(data,rr); l$report=rep;
#	n = unlist(strsplit(file, "/"))
#	name = paste(odir, "P", substr(n[length(n)], 0, nchar(n[length(n)])-4), "_FinalReport_ADM2.txt", sep="")
#	p=l$report[l$report[,4]>0,]; if(length(p)==8) { p = rbind(p,p); }
#	write.table(p, file=name, sep="\t", row.names=F, col.names=F, quote=F)
#	name = paste(odir, "N", substr(n[length(n)], 0, nchar(n[length(n)])-4), "_FinalReport_ADM2.txt", sep="")
#	n=l$report[l$report[,4]<=0,]; if(length(n)==8) { n = rbind(n,n); }
#	write.table(n, file=name, sep="\t", row.names=F, col.names=F, quote=F)
#	return(l)
#}
