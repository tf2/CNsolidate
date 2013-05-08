`MakeData`<-function(TotalTable, infoPos.StartEnd)
{
	if (infoPos.StartEnd) {
	MetaTable<-TotalTable[,1:4]
	} else
	{
	MetaTable<-TotalTable[,1:3]
	}
	if (infoPos.StartEnd) {
	NumericTable<-TotalTable[,5:dim(TotalTable)[2]]
	} else
	{
	NumericTable<-TotalTable[,4:dim(TotalTable)[2]]
	}
	NExp<-(dim(NumericTable)[2])/2
	TableSeg<-as.matrix(NumericTable)
	SummaryData<-c()

	#for (i in 1:NExp)
	#{
	i = 1
		segdata<-as.numeric(TableSeg[,i])
            NData<-NumericTable
		startseg<-c(1,(1+which(diff(segdata)!=0)))
		endseg<-c(which(diff(segdata)!=0),length(segdata))
            sdvec<-c()
            for (j in 1:length(startseg))
            {
                 sdvec<-c(sdvec,sd(NData[startseg[j]:endseg[j]]))
            }
		SummaryData<-rbind(SummaryData,cbind(rep.int(i,length(startseg)),startseg,endseg,segdata[startseg],sdvec))
	#}

	DataList<-list()
	DataList$SummaryData<-SummaryData
	DataList$MetaTable<-MetaTable
	DataList
}


