`GADA_s` <-
function(data, file, idir, odir, alf = 0.1, t = 3, minL = 5, minR=0.2) {
		
		my = cbind(data[,1:3], data[,4])
		my = my[complete.cases(my),]
		my = my[my[,1]!=0,]
		my = my[my[,1]<23,]
		my = my[order(my[,1],my[,2], my[,3]),]
		arrayData=cbind(my[,1:2], my[,4])

		GenInfo = cbind(1:length(arrayData[,1]), arrayData[,1], arrayData[,2])
		arrayData = arrayData[,3]
		colnames(GenInfo) = c("probe", "chr", "pos")
		GenInfo = as.data.frame(GenInfo)
		
		dataSim<-setupGADAgeneral(arrayData,gen.info=GenInfo)
		step1<-SBL(dataSim, sigma2=dLRs(arrayData), aAlpha=alf)
		step2<-BackwardElimination(step1,T=t,MinSegLen=minL)
		res = summary(step2,print=FALSE)
		res = res[complete.cases(res),]
		res = res[res[,6]!=0,]
		res = res[abs(res[,4])>minR,]
		for(x in 1:length(res[,1])) {
			t = data[data[,1] == res[x,5] & data[,2]==res[x,2],3]
			res[x,2] = t[1]
		}
		
		# GADA is strange - it does some scaling - not median substraction (something else)
		#u = unique(data[,1])
		#for(x in 1:length(u)) {
		#	m = median(data[data[,1]==u[x],4])
		#	
		#	res[as.integer(res[,5])==u[x],4] = res[as.integer(res[,5])==u[x],4] +m
		#}

		n = unlist(strsplit(file, "/"))
		name = paste(odir, substr(n[length(n)], 0, nchar(n[length(n)])-4), "_FinalReport_GADA.txt", sep="")
		write.table(cbind(as.integer(res[,5]), as.integer(res[,1]), as.integer(res[,2]), res[,4], res[,3]), file=name, sep="\t", quote=F, col.names=F, row.names=F)
		
return(res)
}

