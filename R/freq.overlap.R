`freq.overlap` <- function(set=NULL, s2=1) {
	set2 = NULL; r = NULL; n = unlist(strsplit(set$files$file, "/"));
	report = paste(set$files$odir, "/", substr(n[length(n)], 0, nchar(n[length(n)])-4), "_FinalReport.txt", sep="")
	set1 = read.table(report);
	if(is.null(s2)) { data(cnv2.1_separate); set2=cnv2.1_separate; }
	if(s2==1) { data(cnv2.1_separate); set2=cnv2.1_separate; }
	if(s2==2) { data(cnv2.1); set2=cnv2.1; } 
	si1 = length(set1[,1]); si2 = length(set2[,1]); 
	set1=set1[set1[,2]<=set1[,3],]; set2=set2[set2[,2]<=set2[,3],];
	if( si1>length(set1[,1]) | si2>length(set2[,1]) ) { warning("i removed some rubbish rows - where start was larger then stop!") }
	if(length(set1)>0 & length(set2)>0) {
	set1 = set1[order(set1[,1], set1[,2],set1[,3]),];
	set2 = set2[order(set2[,1], set2[,2],set2[,3]),];
	u = unique(set1[,1])
	for(x in 1:length(u)) {
		s1 = set1[as.character(set1[,1])==u[x],]
		s2 = set2[as.character(set2[,1])==u[x],]
		r1 = vector(length=length(s1[,1]));
		r2 = vector(length=length(s1[,1]));
		rf1 = vector(length=length(s1[,1]));
		rf2 = vector(length=length(s1[,1]));
		rf3 = vector(length=length(s1[,1]));
		rty = vector(length=length(s1[,1]));
		if(length(s1)>0 & length(s2)>0) {
		if(is.null(dim(s1))) { s1 = rbind(s1,s1); }
		if(is.null(dim(s2))) { s2 = rbind(s2,s2); }
			res <- .C("mapfrequency",
				"a" = as.integer(s1[,2]) 
				,"b" = as.integer(s1[,3])
				,"c" = as.integer(s2[,2]) 
				,"d" = as.integer(s2[,3])
				,"fre1" = as.double(s2[,5]) 
				,"fre2" = as.double(s2[,8]) 
				,"fre3" = as.double(s2[,11]) 
				,"r1" = as.double(r1) 
				,"r2" = as.double(r2)
				,"rf1" = as.double(rf1) 
				,"rf2" = as.double(rf2) 
				,"rf3" = as.double(rf3)
				,"ty" = as.double(s2[,13]) 
				,"rty" = as.double(rty)
				,"size1" = as.integer(length(s1[,1]))
				,"size2" = as.integer(length(s2[,2])) 
				,"PACKAGE" = "CNsolidate")
			r = rbind(r, cbind(res$r1, res$r2, res$rf1, res$rf2, res$rf3, res$rty))
		} else {
			r = rbind(r, cbind(r1, r2, -1, -1, -1, -1))
		}
	}
	rset = cbind(set1,r)
	write.table(rset, file=report, sep="\t", row.names=F, col.names=F, quote=F)
	invisible(rset)
	} else {
		 warning("everything was rubbish rows - where start was larger then stop!")
	}
}
