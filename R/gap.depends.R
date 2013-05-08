`gap.depends` <- function(set=NULL) {
	n = unlist(strsplit(set$files$file, "/"))
    	rep = paste(set$files$odir, "/", substr(n[length(n)], 0, nchar(n[length(n)]) - 4), "_FinalReport.txt", sep = "")
	print(rep)
	r = read.table(rep); r = r[order(r[,1], r[,2], r[,3]),]
	d = read.table(set$files$file);  d = d[order(d[,1], d[,2], d[,3]),]
	cent = paste(system.file("extdata", package="CNsolidate"), "/centromeres_hg19.txt", sep="")
	cens = read.table(cent, sep="\t");  cens = cens[order(cens[,1], cens[,2], cens[,3]),]
	
	ov = uber.overlap(r, cens)
	ind = length(ov[1,])-3
	fine = ov[ov[,ind]==0,]
	notfine = ov[ov[,ind]>0,]
	
	fineagain = NULL
	if(length(notfine[, 1])>0) {
		for(x in 1:length(notfine[,1])) {
			chr = notfine[x,1]
			start = notfine[x,2]
			stop = notfine[x,3]
			dd = d[d[,1]==chr,]
			cenrow = cens[as.numeric(as.character(notfine[x, length(notfine[x,])])),]
			newstop = cenrow[1,2]
			ddd = dd[dd[,3]<newstop,]
			newstop = ddd[which.min(abs(ddd[,3]-newstop)),3]
		
			if(length(newstop)==0) {
				newstop = cenrow[1,3]
				ddd = dd[dd[,3]<newstop,]
				newstop = ddd[which.min(abs(ddd[,3]-newstop)),3]
			}
			
			ind1 = notfine[x,6]
			n_ind2 = which.min(abs(ddd[,3]-newstop))
			m1 = mean(dd[ind1:n_ind2,4])
		
			newstart = cenrow[1,3]
			ddd = dd[dd[,2]>newstart,]
			newstart = ddd[which.min(abs(ddd[,2]-newstart)),2]

			if(length(newstart)==0) {
				newstart = cenrow[1,2]
				ddd = dd[dd[,2]>newstart,]
				newstart = ddd[which.min(abs(ddd[,2]-newstart)),2]
			}
			
			ind2 = notfine[x,7]
			n_ind1 = which.min(abs(dd[,2]-newstart))
			m2 = mean(dd[n_ind1:ind2,4])
		
			fineagain=rbind(fineagain, cbind(chr, start, newstop, m1, (n_ind2-ind1)+1, ind1, n_ind2))
			fineagain=rbind(fineagain, cbind(chr, newstart, stop, m2, (ind2-n_ind1)+1, n_ind1, ind2))
		}
		en = ind-3
		colnames(fineagain) = colnames(fine[,1:en])
	}

	en = ind-3
	r = rbind(fine[,1:en], fineagain)
	r = r[r[, 5] >= set$combine.settings$minProbe & abs(r[,4]) >= set$combine.settings$absRatio, ]
	r=r[order(r[,1], r[,2], r[,3]),]
	
	write.table(r, file=rep, sep="\t", row.names=F, col.names=F, quote=F)
	
invisible(r)	
}
