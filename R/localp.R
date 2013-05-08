`localp` <- function(set=NULL) {
	n = unlist(strsplit(set$files$file, "/")); tsig=vector(); vals=vector();
	report = paste(set$files$odir, "/", substr(n[length(n)], 0, nchar(n[length(n)])-4), "_FinalReport.txt", sep="");
	data=read.table(set$files$file); rep=read.table(report);
	pin = 1; l1=list(); pin = 1; u=unique(rep[,1]);
	for(y in 1:length(u)) {
	d=data[data[,1]==u[y],]; v=rep[rep[,1]==u[y],]; dat=d;
		for(x in 1:length(v[,1])) {
			t=v[x,]; n=t[5]; ind1=which(as.numeric(dat[,2])==as.numeric(t[2])); ind2=which(as.numeric(dat[,3])==as.numeric(t[3]));
			if(length(ind1)>0 & length(ind2)>0) {
				i1 = as.numeric(ind1[1]); i2 = as.numeric(ind2[length(ind2)]); dat=dat[-(i1:i2),];
			}
		}
		su = 0;
		for(x in 1:length(v[,1])) {
			t=v[x,];  n=t[5]; ind1=t[6]; ind2=t[7];
			p1= d[as.numeric(ind1[1]):as.numeric(ind2[length(ind2)]),4];
			ind1=which(as.numeric(d[,2])==as.numeric(t[2])); ind2=which(as.numeric(d[,3])==as.numeric(t[3]));
			ind1 = ind1-su; ind2 = ind1+(n*set$combine.settings$tfac)-1; #ind1 = (ind1[1]-(n*set$combine.settings$tfac)); 
			su=su+n;
			ind1=ind1[1]; ind2=ind2[length(ind2)];
			if(ind1<1) { 
				tem = abs(ind1); ind2=ind2+tem; 
				if (ind2 > length(dat[, 1])) {
                	ind2 = length(dat[, 1]);
                } 
                ind1=1;
            }
			if(ind2>length(dat[,1])) { 
				temp=ind2-length(dat[,1]); ind1=ind1-temp; 
				 if (ind1 < 1) {
                 	ind1=1;
                 }
				ind2=length(dat[,1]); 
			}
			p2=dat[as.numeric(ind1[1]):as.numeric(ind2[length(ind2)]),4];
			val = t.test(p1, p2)$p.value; ll = list();
			ll$value = val; ll$t=t; ll$p1 = p1; ll$p2 = p2;
			l1=c(l1,ll); tsig[pin]=val; vals[pin]=paste(p1,collapse=","); pin=pin+1; 
		}
	}
	#write.table(cbind(rep,tsig, vals), file=report, sep="\t", row.names=F, col.names=F, quote=F);
	write.table(cbind(rep,tsig), file=report, sep="\t", row.names=F, col.names=F, quote=F);
	invisible(l1)
}	


#`localp` <- function(file, odir, fac=2) {
#	n = unlist(strsplit(file, "/")); tsig=vector(); vals=vector();
#	report = paste(odir, substr(n[length(n)], 0, nchar(n[length(n)])-4), "_FinalReport.txt", sep="");
#	data=read.table(file); rep=read.table(report);
#	pin = 1; l1=list(); pin = 1; u=unique(rep[,1]);
#	for(y in 1:length(u)) {
#	d=data[data[,1]==u[y],]; v=rep[rep[,1]==u[y],]; dat=d;
#		for(x in 1:length(v[,1])) {
#			t=v[x,]; n=t[5]; ind1=which(as.numeric(dat[,2])==as.numeric(t[2])); ind2=which(as.numeric(dat[,3])==as.numeric(t[3]));
#			if(length(ind1)>0 & length(ind2)>0) {
#				i1 = as.numeric(ind1[1]); i2 = as.numeric(ind2[length(ind2)]); dat=dat[-(i1:i2),];
#			}
#		}
#		su = 0;
#		for(x in 1:length(v[,1])) {
#			t=v[x,];  n=t[5]; ind1=t[6]; ind2=t[7];
#			p1= d[as.numeric(ind1[1]):as.numeric(ind2[length(ind2)]),4];
#			ind1=which(as.numeric(d[,2])==as.numeric(t[2])); ind2=which(as.numeric(d[,3])==as.numeric(t[3]));
#			ind1 = ind1-su; ind2 = ind1+(n*fac)-1; #ind1 = (ind1[1]-(n*fac)); 
#			su=su+n;
#			ind1=ind1[1]; ind2=ind2[length(ind2)];
#			if(ind1<1) { 
#				tem = abs(ind1); ind2=ind2+tem; 
#				if (ind2 > length(dat[, 1])) {
#               	ind2 = length(dat[, 1]);
#                } 
#               ind1=1;
#            }
#			if(ind2>length(dat[,1])) { 
#				temp=ind2-length(dat[,1]); ind1=ind1-temp; 
#				 if (ind1 < 1) {
#                	ind1=1;
#                 }
#				ind2=length(dat[,1]); 
#			}
#			p2=dat[as.numeric(ind1[1]):as.numeric(ind2[length(ind2)]),4];
#			val = t.test(p1, p2)$p.value; ll = list();
#			ll$value = val; ll$t=t; ll$p1 = p1; ll$p2 = p2;
#			l1=c(l1,ll); tsig[pin]=val; vals[pin]=paste(p1,collapse=","); pin=pin+1; 
#		}
#	}
#	write.table(cbind(rep,tsig, vals), file=report, sep="\t", row.names=F, col.names=F, quote=F);
#	return(l1)
#}	
