`FastCall`<-function(MatrixTot, Factor = 2, thr0 = 0.10, Nclass = 4, sample.names = NULL)
{
	time.Started <- proc.time()
	Nclones<-dim(MatrixTot)[1]
	ix.tmp<-round(runif(100)*(Nclones-1))+1
	if (abs(mean(as.numeric(MatrixTot[ix.tmp,4]), na.rm=TRUE)) < 2) {
		infoPos.StartEnd = FALSE
		nexp <- (length(MatrixTot[1,])-3)/2
	} else
	{
		infoPos.StartEnd = TRUE
		nexp <- (length(MatrixTot[1,])-4)/2
	}
	if (is.null(sample.names)) {sample.names <- seq.int(1, nexp)}
	Call.List<-MakeData(MatrixTot, infoPos.StartEnd)
	MetaData<-Call.List$MetaTable
	SummaryData<-Call.List$SummaryData
	mdata<-SummaryData[,4]
	ResultsEM<-EMFastCall(mdata,thr0,Factor)
	muvec<-ResultsEM$muvec
	sdvec<-ResultsEM$sdvec
	prior<-ResultsEM$prior
	iter<-ResultsEM$iter
	P0<-PosteriorP(mdata,muvec,sdvec,prior)
	out <- LabelAss(P0, mdata, Nclass)
	FinalTable<-DataRecomp(SummaryData,MetaData,out)
	timeFC <- (proc.time() - time.Started)[1]
	
	
	
	C.Call <- paste("Call sample ", sample.names)
	C.Prob <- paste("Prob sample", sample.names)
	C.Log <- paste("Logratio sample", sample.names)
	C.Seg <- paste("Seg sample", sample.names)
	if (nexp == 1) {
		C.Call <- "Call"
		C.Prob <- "Prob"
		C.Log <- "Logratio"
		C.Seg <- "Seg"
	}
	
	if (infoPos.StartEnd) {
		colnames(FinalTable) <- c("ProbeID","Chromosome","Start Position","End Position",C.Call,C.Prob)
		MatrixSeg<-MatrixTot
		colnames(MatrixSeg) <- c("ProbeID","Chromosome","Start Position","End Position",C.Log,C.Seg)
	} else
	{
		colnames(FinalTable) <- c("ProbeID","Chromosome","Position",C.Call,C.Prob)
		MatrixSeg<-MatrixTot
		colnames(MatrixSeg) <- c("ProbeID","Chromosome","Position",C.Log,C.Seg)
	}
	
	Results<-list()
	Results$Call<-FinalTable
	Results$Seg<-MatrixSeg
	Results$timeFC <- timeFC
	
	return(Results)
	
}

