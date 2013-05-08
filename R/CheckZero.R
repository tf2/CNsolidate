`CheckZero`<-function(normaldataMat,mdata,lvec,uvec)
{
Ncount<-dim(normaldataMat)[2]
tmp<-rowSums(normaldataMat)
indtmp<-which(tmp==0)
if (length(indtmp!=0))
{
	for (k in 1:length(indtmp))
	{
		dataMattmp<-rep(0,Ncount)
		mtmp<-mdata[indtmp[k]]
		for (ki in 1:Ncount)
			{
				l<-lvec[ki]
				u<-uvec[ki]
				dataMattmp[ki]<-(mtmp<=u)*(mtmp>=l)
			}
		normaldataMat[indtmp[k],]<-dataMattmp
	}
}
return(normaldataMat)
}