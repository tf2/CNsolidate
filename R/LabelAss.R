`LabelAss`<-function(P0, mdata, Nclass)
{
	CallResults<-rep(0,length(mdata))
	ProbResults<-rep(0,length(mdata))
	
#indcall<-max.col(P0)
	if (Nclass == 4) 
	{
		P0.Class <- cbind((P0[,1]+P0[,2]), P0[,3:5])
		indcall<-max.col(P0.Class)
	CallResults[indcall==1]<- -1
	CallResults[indcall==2] <- 0
	CallResults[indcall==3]<- 1
	CallResults[indcall==4]<- 2
#	CallResults[indcall==5]<- 2
	for (i in 1:length(mdata))
	{
	ProbResults[i]<-P0.Class[i,indcall[i]]
	}
	}
	if (Nclass == 3) 
	{
		P0.Class <- cbind((P0[,1]+P0[,2]), P0[,3], (P0[,4]+P0[,5]))
		indcall<-max.col(P0.Class)
		CallResults[indcall==1]<- -1
		CallResults[indcall==2] <- 0
		CallResults[indcall==3]<- 1
#		CallResults[indcall==4]<- 2
#	CallResults[indcall==5]<- 2
		for (i in 1:length(mdata))
		{
			ProbResults[i]<-P0.Class[i,indcall[i]]
		}
	}
	
Results<-c()
Results<-cbind(CallResults,ProbResults)
Results
}


