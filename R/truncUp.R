`truncUp`<-function(muvec,sdvec,thr0,Factor)
{
	lvec<-c()
	uvec<-c()
	sigmaNeut<-Factor*sdvec[3]
	
	if (sigmaNeut<thr0)
	{
		sigmaNeut<-thr0
	}
	sigmaGain<-Factor*sdvec[4]
	sigmaLoss<-Factor*sdvec[2]
	uvec4<-muvec[4]+sigmaGain
	uvec1<-muvec[2]-sigmaLoss
	if (uvec4>0.9)
	{
	uvec4<-0.9
	}
	if (uvec1< -1.1)
	{
	uvec1<- -1.1
	}
	lvec[1]<- -16
	uvec[1]<- uvec1
	lvec[2]<- uvec1
	uvec[2]<- -sigmaNeut
	lvec[3]<- -sigmaNeut
	uvec[3]<-  sigmaNeut
	lvec[4]<-  sigmaNeut
	uvec[4]<-  uvec4
	lvec[5]<-  uvec4
	uvec[5]<- 16
	
Result<-list()
Result$l<-lvec
Result$u<-uvec
Result
}


