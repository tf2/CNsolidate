`meeting.c.o` <- function(file, idir, odir, conNum=1, methods = c("safe", "walk", "SMUG", "GADA", "SMAP", "dna.copy")) {

	command = "cat";

	fName = paste(odir, substr(file, 0, nchar(file)-4), "_Alltemp.txt", sep="")
	finName = paste(odir, substr(file, 0, nchar(file)-4), "_Allfeatures.txt", sep="")
	finRep = paste(odir, substr(file, 0, nchar(file)-4), "_FinalReport.txt", sep="")

	for(x in 1:length(methods)) {
		name1 = paste(odir, "N", substr(file, 0, nchar(file)-4), "_FinalReport_", methods[x], ".txt", sep="")
		name2 = paste(odir, "P", substr(file, 0, nchar(file)-4), "_FinalReport_", methods[x], ".txt", sep="")
		d1 = read.table(name1, sep="\t")
		d2 = read.table(name2, sep="\t")
		d1=cbind(d1[,1:4], as.character(methods[x]))
		d2=cbind(d2[,1:4], as.character(methods[x]))
		write.table(d1, file=name1, sep="\t", row.names=F, col.names=F, quote=F)
		write.table(d2, file=name2, sep="\t", row.names=F, col.names=F, quote=F)
		command = paste(command, name1, name2, sep=" ")
	}
	
	command = paste(command, ">", finName, sep=" ")
	system(command)
	#command = paste("cut -f 1-4 ", fName, ">", finName, sep="")
	#system(command)
	
	d = read.table(finName, sep="\t")
	d[,1] = as.numeric(d[,1])
	d=d[order(d[,1], d[,2]),]
	write.table(d, file=finName, sep="\t", row.names=F, col.names=F, quote=F)
	
	n1 = paste(odir, "N", file, sep="")
	n2 = paste(odir, "P", file, sep="")
	tfile = paste(odir, file, sep="")
	command = paste("cat", n1, n2, ">", tfile, sep=" ")
	system(command)
	
	command = paste("rm", n1, n2, sep=" ")
	system(command)
	
	d = read.table(tfile, sep="\t")
	d[,1] = substr(d[,1], 4,10)
	d[,1] = as.numeric(d[,1])
	d=d[order(d[,1], d[,4]),]
	write.table(d, file=tfile, sep="\t", row.names=F, col.names=F, quote=F)
	
	Jpath =  system.file("java", package="CNsolidate")
	setwd(Jpath)
	command = paste("java -Xmx1600m makeObject -f ", tfile, " -r ", finName, " -n ", conNum, " >", fName, sep="")
	system(command)
	
	command = paste("java -Xmx1600m Fmerge -i ", tfile, " -o ", fName, " >", finRep, sep="")
	system(command)
	
	command = paste("java -Xmx1600m finishOff -r ", finRep, sep="")
	system(command)
}