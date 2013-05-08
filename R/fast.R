`fast` <- function(set=NULL) {
	data = set$data
	if(is.null(data)) { data = read.table(set$files$file); }
	fname = paste(substr(set$files$file,0,nchar(set$files$file)-4), "_FinalReport_fast.txt", sep="")
 	fa <- function(data, set=NULL) {
 		m = cbind(1:length(data[,1]), data[,1], data[,2], data[,4], runmed(data[,4],301))
 		uni = unique(m[,2]); A = list(); l = NULL
 			for(x in 1:length(uni)) {
				if(uni[x] < 23) {
					d = m[m[,2]==uni[x],]
					d[,4:5] = abs(d[,4:5])
					dd = data[data[,1]==uni[x],]
					print(uni[x])
					A[[x]] = cbind(dd[,1], dd[,2], dd[,3], dd[,4], FastCall(d)$Call[,4:5])
				}
 			}	
			for(x in 1:length(A)) {
				l = rbind(l, A[[x]])	
			}
		name = paste(substr(set$files$file,0,nchar(set$files$file)-4), "_fast.temp", sep="")
		fname = paste(substr(set$files$file,0,nchar(set$files$file)-4), "_FinalReport_fast.txt", sep="")
		ll= data.frame(as.integer(l[,1]), as.integer(l[,2]), as.integer(l[,3]), l[,-(1:3)])
		write.table(ll, file=name, sep="\t", row.names=F, col.names=F, quote=F)
		Jpath =  system.file("java", package="CNsolidate")
		setwd(Jpath); command = paste("java -Xmx1600m fiFast -f ", name, " -o ", fname, sep=""); system(command);
		command = paste("rm ", name,sep=""); system(command);
	}
	fa(data[data[,4]>0,], set); f1 = read.table(fname);
	fa(data[data[,4]<0,], set); f2 = read.table(fname);
	f = rbind(f1,f2); f = f[order(f[,1],f[,2],f[,3]),];
	write.table(f, file=fname, sep="\t", row.names=F, col.names=F, quote=F)
	invisible(f)
}



#`rFASTCALL` <- function(data, file, idir, odir) {
#
#	m = cbind(1:length(data[,1]), data[,1], data[,2], data[,4], runmed(data[,4],301))
# 	name = paste(substr(file,0,nchar(file)-4), "_fast.temp", sep="")#
#	fname = paste(substr(file,0,nchar(file)-4), "_FinalReport_fast.txt", sep="")
# 	uni = unique(m[,2])
# 	A = list()
# 	for(x in 1:length(uni)) {
#		if(uni[x] < 23) {
#			d = m[m[,2]==uni[x],]
#			d[,4:5] = abs(d[,4:5])
#			dd = data[data[,1]==uni[x],]
#			print(uni[x])
#			A[[x]] = cbind(dd[,1], dd[,2], dd[,3], dd[,4], FastCall(d)$Call[,4:5])
#		}
# 	}
# 	l = NULL	
#	for(x in 1:length(A)) {
#		l = rbind(l, A[[x]])	
#	}
#	l[,1] = as.integer(l[,1])
#	l[,2] = as.integer(l[,2])
#	l[,3] = as.integer(l[,3])
#	write.table(l, file=name, sep="\t", row.names=F, col.names=F, quote=F)
#	Jpath =  system.file("java", package="CNsolidate")
#	setwd(Jpath)
#	command = paste("java -Xmx1600m fiFast -f ", name, " -o ", fname, sep="")
#	system(command)
#	command = paste("rm ", name,sep="")
#	system(command)
#	return(l)
#}