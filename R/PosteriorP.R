`PosteriorP`<-function(mdata,muvec,sdvec,prior)
{
		Ncount<-length(muvec)
		Posterior <- matrix(0,nrow=length(mdata),ncol=length(muvec))
		deno <- 0
		normaldataMat<-c()

		for (j in 1:Ncount)
		{
			moy<-muvec[j]
			sdev<-sdvec[j]
			normaldataMat<-cbind(normaldataMat,dnorm(mdata,mean=moy,sd=sdev))
		}
		tauMat<-c()
		for (j in 1:Ncount)
		{
			tauMat<-cbind(tauMat,c(prior[j]*(t(normaldataMat[,j]))))		}
			
		deno <- rowSums(tauMat)
		Posterior<- tauMat/deno
		Posterior
}
