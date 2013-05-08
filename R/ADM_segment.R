`ADM_segment` <- function(dd, t=5, rd=2) {
	dd[dd[,4]==0,4]=0.0001;
	rr1=dd; rr1[rr1[,4]>0,6]=0; rr1[rr1[,4]>0,7]=max(rr1[rr1[,4]>0,7]);
	rrr <- .C("_EADM"
			,"d" = as.double(rr1[,6])
			,"dat" = as.double(rr1[,6])
			,"prob" = as.double(rr1[,7])
		 	,"size" = as.integer(length(rr1[,1]))
		 	,"t" = as.double(t)
		 	,"rd" = as.double(rd)
		 	,"PACKAGE" = "CNsolidate")
	rrrr=cbind(rr1[,1], rr1[,2], rr1[,3], rr1[,4], rr1[,5], rrr$dat, rrr$prob, rrr$d); rr1=rrrr[rrrr[,8]!=0,];
	rr2=dd; rr2[rr2[,4]<0,6]=0; rr2[rr2[,4]<0,7]=max(rr2[rr2[,4]<0,7]);
	rrr <- .C("_EADM"
			,"d" = as.double(rr2[,6])
			,"dat" = as.double(rr2[,6])
			,"prob" = as.double(rr2[,7])
		 	,"size" = as.integer(length(rr2[,1]))
		 	,"t" = as.double(t)
		 	,"rd" = as.double(rd)
		 	,"PACKAGE" = "CNsolidate")
	rrrr=cbind(rr2[,1], rr2[,2], rr2[,3], rr2[,4], rr2[,5], rrr$dat, rrr$prob, rrr$d); rr2=rrrr[rrrr[,8]!=0,]; r=rbind(rr1,rr2);
	return(r[order(r[,1],r[,2],r[,3]),])
}