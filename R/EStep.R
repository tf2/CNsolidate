`EStep`<-function(mdata,muvec,sdvec,prior,lvec,uvec)
{
	Ncount<-length(muvec)
	taux <- matrix(0,nrow=length(mdata),ncol=Ncount)
	deno <- 0
	normaldataMat<-c()
	for (g in 1:Ncount)
	{
		mu<-muvec[g]
		sdev<-sdvec[g]
		l<-lvec[g]
		u<-uvec[g]
		if (sum((mdata<=u)*(mdata>=l))!=0)
		{
		normaldataMat<-cbind(normaldataMat,c(gfct(t(mdata),mu,sdev,l,u)))
		} else
		{
		normaldataMat<-cbind(normaldataMat,rep(0,length(mdata)))
		}
	}
	normaldataMat<-CheckZero(normaldataMat,mdata,lvec,uvec)
	tauMat<-c()
	for (j in 1:Ncount)
	{
		tauMat<-cbind(tauMat,c(prior[j]*(t(normaldataMat[,j]))))
	}
	deno <- rowSums(tauMat)
	taux<- tauMat/deno
taux
}
