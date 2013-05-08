`GADA` <- function(set=NULL) {
	data = set$data
	if(is.null(data)) { data = read.table(set$files$file); }
	dat = cbind(data[,1:3], data[,4]); dat = dat[complete.cases(dat),]
	dat = dat[dat[,1]!=0,]; dat = dat[dat[,1]<23,];
	dat = dat[order(dat[,1],dat[,2], dat[,3]),]
	arrayData=cbind(dat[,1:2], dat[,4])

	GenInfo = cbind(1:length(arrayData[,1]), arrayData[,1], arrayData[,2]); arrayData = arrayData[,3];
	colnames(GenInfo) = c("probe", "chr", "pos"); GenInfo = as.data.frame(GenInfo);
		
	dataSim<-setupGADAgeneral(arrayData,gen.info=GenInfo)
	step1<-SBL(dataSim, sigma2=dLRs(arrayData), aAlpha=set$GADA$alpha)
	step2<-BackwardElimination(step1,T=set$GADA$T,MinSegLen=set$GADA$mL)
	res = summary(step2,print=FALSE); res = res[complete.cases(res),];

		for(x in 1:length(res[,1])) {
			t = data[data[,1] == res[x,5] & data[,2]==res[x,2],3]
			res[x,2] = t[1]
		}
	
	res = res[res[,6]!=0,]
	res = res[abs(res[,4])>set$GADA$mR,]
	n = unlist(strsplit(set$files$file, "/"))
	dat = data.frame(as.integer(res[,5]), as.integer(res[,1]), as.integer(res[,2]), res[,4], res[,3])
	
	#name = paste(set$files$odir, "P", substr(n[length(n)], 0, nchar(n[length(n)])-4), "_FinalReport_GADA.txt", sep="")
	#write.table(dat[dat[,4]>0,], file=name, sep="\t", quote=F, col.names=F, row.names=F)
	
	name = paste(set$files$odir, "/", substr(n[length(n)], 0, nchar(n[length(n)])-4), "_FinalReport_GADA.txt", sep="")
	write.table(dat, file=name, sep="\t", quote=F, col.names=F, row.names=F)	

invisible(res)
}


#`GADA` <- function(data, file, idir, odir, alf = 0.1, t = 3, minL = 5, minR=0.2) {
#		
#		my = cbind(data[,1:3], data[,4])
#		my = my[complete.cases(my),]
#		my = my[my[,1]!=0,]
#		my = my[my[,1]<23,]
#		my = my[order(my[,1],my[,2], my[,3]),]
#		arrayData=cbind(my[,1:2], my[,4])
#
#		GenInfo = cbind(1:length(arrayData[,1]), arrayData[,1], arrayData[,2])
#		arrayData = arrayData[,3]
#		colnames(GenInfo) = c("probe", "chr", "pos")
#		GenInfo = as.data.frame(GenInfo)
#		
#		dataSim<-setupGADAgeneral(arrayData,gen.info=GenInfo)
#		step1<-SBL(dataSim, sigma2=dLRs(arrayData), aAlpha=alf)
#		step2<-BackwardElimination(step1,T=t,MinSegLen=minL)
#		res = summary(step2,print=FALSE)
#		res = res[complete.cases(res),]
#
#		for(x in 1:length(res[,1])) {
#			t = data[data[,1] == res[x,5] & data[,2]==res[x,2],3]
#			res[x,2] = t[1]
#		}
#	
#	    res = res[res[,6]!=0,]
#		res = res[abs(res[,4])>minR,]
#		n = unlist(strsplit(file, "/"))
#		dat = data.frame(as.integer(res[,5]), as.integer(res[,1]), as.integer(res[,2]), res[,4], res[,3])
#		name = paste(odir, "P", substr(n[length(n)], 0, nchar(n[length(n)])-4), "_FinalReport_GADA.txt", sep="")
#		write.table(dat[dat[,4]>0,], file=name, sep="\t", quote=F, col.names=F, row.names=F)
#		name = paste(odir, "N", substr(n[length(n)], 0, nchar(n[length(n)])-4), "_FinalReport_GADA.txt", sep="")
#		write.table(dat[dat[,4]<0,], file=name, sep="\t", quote=F, col.names=F, row.names=F)
#		
#return(res)
#}
