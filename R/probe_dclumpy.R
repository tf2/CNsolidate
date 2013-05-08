probe_dclumpy <- function(file, odir, mlen = 3, abCut = 0.4, sizeCut = 5000) {
	n = unlist(strsplit(file, "/"))
	report = paste(odir, substr(n[length(n)], 0, nchar(n[length(n)])-4), "_FinalReport.txt", sep="")
	rep = read.table(report, sep="\t"); nsc = vector(); data = read.table(file);		
	u = unique(rep[,1]); v = vector(); vv=vector(); count = 1; re = rep
	for(x in 1:length(u)) {	
		d = data[data[,1]==u[x],]; r = rep[rep[,1]==u[x],];
		for(y in 1:length(r[,1])) {
			row = r[y,];
			sts = d[as.numeric(row[6]):as.numeric(row[7]),]
			s = max(diff(sts[,2]))/(max(sts[length(sts[,1]),2]-sts[1,2]))
			vv[count] = 0		
			if(s=="NaN") { s=0; }
			if(s > 0.85) {
				di = diff(sts[,2]); ind = which.max(di); size = max(di);
				si = c(abs(1),abs(length(sts[,1]))); side= which.min(c(abs(1-ind),abs(length(sts[,1])-ind)));
				val = abs(si[side]-ind)
				if (side==1) { val=val+1; }
				num = 0; me = row[4]/2;
				if(me < abCut) {
					if (me < 0) { me = -abCut; }
					if (me >= 0) { me = abCut; }	
				}
				if(side==1) {	
					for(z in 1:val) {
						if(me < 0 & sts[z,4]< me) { num=num+1; }
						if(me >=0 & sts[z,4]> me) { num=num+1; }	
					}
				}
				if(side==2) {	
					for(z in length(sts[,1]):val) {
						if(me < 0 & sts[z,4]< me) { num=num+1; }
						if(me >=0 & sts[z,4]> me) { num=num+1; }	
					}
				}
				if(num < mlen & size > sizeCut) {	
					vv[count] = 1;
					if(side == 1) {
						re[count,6] = re[count,6]+val;
						re[count,2] = d[as.numeric(re[count,6]),2];
						re[count,3] = d[as.numeric(re[count,7]),3];
						re[count,4] = mean(d[as.numeric(re[count,6]):as.numeric(re[count,7]),4]);
						re[count,5] = length(d[as.numeric(re[count,6]):as.numeric(re[count,7]),4]);
					}	
					if(side == 2) {
						re[count,7] = re[count,7]-val;
						re[count,2] = d[as.numeric(re[count,6]),2];
						re[count,3] = d[as.numeric(re[count,7]),3];
						re[count,4] = mean(d[as.numeric(re[count,6]):as.numeric(re[count,7]),4]);
						re[count,5] = length(d[as.numeric(re[count,6]):as.numeric(re[count,7]),4]);
					}
					if((re[count,7]-re[count,6])<3) { vv[count] = 2 }
				}
			}
		v[count] = s; count=count+1;
		}
	}
	return(re)
}	