`rFASTCALL` <- function(data, file, idir, odir) {

 	m = cbind(1:length(data[,1]), data[,1], data[,2], data[,4], runmed(data[,4],301))
 	name = paste(substr(file,0,nchar(file)-4), "_fast.temp", sep="")
	fname = paste(substr(file,0,nchar(file)-4), "_FinalReport_fast.txt", sep="")
 	uni = unique(m[,2])
 	A = list()
 	for(x in 1:length(uni)) {
		if(uni[x] < 23) {
			d = m[m[,2]==uni[x],]
			d[,4:5] = abs(d[,4:5])
			dd = data[data[,1]==uni[x],]
			print(uni[x])
			A[[x]] = cbind(dd[,1], dd[,2], dd[,3], dd[,4], FastCall(d)$Call[,4:5])
		}
 	}
 	l = NULL	
	for(x in 1:length(A)) {
		l = rbind(l, A[[x]])	
	}
	l[,1] = as.integer(l[,1])
	l[,2] = as.integer(l[,2])
	l[,3] = as.integer(l[,3])
	write.table(l, file=name, sep="\t", row.names=F, col.names=F, quote=F)
	Jpath =  system.file("java", package="CNsolidate")
	setwd(Jpath)
	command = paste("java -Xmx1600m fiFast -f ", name, " -o ", fname, sep="")
	system(command)
	command = paste("rm ", name,sep="")
	system(command)
	return(l)
}