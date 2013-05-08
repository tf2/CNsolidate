`DataRecomp`<-function(SummaryData,MetaData,out)
{
	ExpClassTotal<-c()
	ExpProbTotal<-c()
	NumExp<-length(unique(SummaryData[,1]))
	for (i in 1:NumExp)
	{
		indExp<-which(SummaryData[,1]==i)
		StartEndMat<-rbind(SummaryData[indExp,2:3])
		OutClass<-out[indExp,1]
		OutProb<-out[indExp,2]
		ExpProb<-c()
		ExpClass<-c()
		for (j in 1:length(indExp))
		{
			ExpProb[StartEndMat[j,1]:StartEndMat[j,2]]<-OutProb[j]
			ExpClass[StartEndMat[j,1]:StartEndMat[j,2]]<-OutClass[j]
		}
	ExpClassTotal<-cbind(ExpClassTotal,ExpClass)
	ExpProbTotal<-cbind(ExpProbTotal,ExpProb)
	}
TotalTableFastCall<-cbind(MetaData,ExpClassTotal,ExpProbTotal)
TotalTableFastCall
}



