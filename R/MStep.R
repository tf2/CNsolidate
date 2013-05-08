`MStep`<-function(mdata,taux,muvec,sdvec)
{
	N <- ncol(taux)
	moy<-c()
	sdev<-c()
	ptmp <- colSums(taux)
	ptmp[ptmp=="NaN"] = 0
	pnew<-c()
	for (j in 1:N)
		{
		if (ptmp[j]!=0)
		{
			moy[j] <- (taux[,j]%*%mdata)/ptmp[j]
			if (sqrt((taux[,j]%*%((mdata-moy[j])^2))/ptmp[j])< 1e-100)
			{ 
				sdev[j]<-sdvec[j]
			} else
			{
			sdev[j]<- sqrt((taux[,j]%*%((mdata-moy[j])^2))/ptmp[j])
			}
			pnew[j]<-ptmp[j]/length(mdata)
		} else
		{
		moy[j] <- muvec[j]
		sdev[j]<- sdvec[j]
		pnew[j]<-1e-06
		}
		}
        pnew<-pnew/sum(pnew)
	if (sdev[1]<sdev[2])
	{
	sdev[1]=sdev[2]
	}
	if (sdev[2]<sdev[3])
	{
	sdev[2]=sdev[3]
	}
	if (sdev[4]<sdev[3])
	{
	sdev[4]=sdev[3]
	}
	if (sdev[5]<sdev[4])
	{
	sdev[5]=sdev[4]
	}
	ParamResult<-list()
	ParamResult$mu<-moy
	ParamResult$sdev<-sdev
	ParamResult$prior<-pnew
ParamResult
}

