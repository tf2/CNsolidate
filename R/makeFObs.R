`makeFObs` <- function(set=NULL, rm=T) {

	file=set$files$file; idir=set$files$idir; odir=set$files$odir;
	methods=NULL; for(x in 1:length(set$algorithms)) { if(set$algorithms[x] == 1) methods=c(methods, names(set$algorithms[x])) }
	conNum=set$combine.settings$Consenues; absCut=set$combine.settings$absRatio; gold = set$files$gold;
	
	command = "cat"; rcommand = "rm";
	
	n = unlist(strsplit(file, "/"))
	finName = paste(odir, "/", substr(n[length(n)], 0, nchar(n[length(n)])-4), "_Allfeatures.txt", sep="")
	finRep = paste(odir, "/", substr(n[length(n)], 0, nchar(n[length(n)])-4), "_FinalReport.txt", sep="")
	for(x in 1:length(methods)) {
		name = paste(odir, "/", substr(n[length(n)], 0, nchar(n[length(n)])-4), "_FinalReport_", methods[x], ".txt", sep="")
		if(length(count.fields(name))>0) {
			d1 = read.table(name, sep="\t")
			if(length(d1[1,])==1) { d1 = t(d1); }
			d1=cbind(as.integer(d1[,1]), as.integer(d1[,2]), as.integer(d1[,3]), as.double(d1[,4]), as.character(methods[x]))
			write.table(d1, file=name, sep="\t", row.names=F, col.names=F, quote=F)
		}
		command = paste(command, name, sep=" "); rcommand = paste(rcommand, name, sep=" ");
	}
	
	command = paste(command, ">", finName, sep=" "); system(command);
	if(rm) {
		system(rcommand)
	}
	
	d = read.table(finName, sep="\t")
	d[,1] = as.numeric(d[,1]); d=d[order(d[,1], d[,2], d[,3]),];
	d=d[abs(d[,4])>absCut,]; d=d[complete.cases(d),];
	write.table(d, file=finName, sep="\t", row.names=F, col.names=F, quote=F)
	
	#command = paste("rm", n1, n2, sep=" ") #system(command)
	
	if(length(count.fields(finName))>0) {
	rep = read.table(finName, sep="\t"); d = read.table(file, sep="\t");
	d=d[order(d[,1], d[,2], d[,3]),]; repdat = NULL; u = unique(rep[,1]); 
	for(x in 1:length(u)) {	
		r=d[d[,1]==u[x],]; rr = rep[rep[,1]==u[x],];
		write.table(r, file=file, sep="\t", row.names=F, col.names=F, quote=F);
		write.table(rr, file=finName, sep="\t", row.names=F, col.names=F, quote=F);
		Jpath =  system.file("java", package="CNsolidate")
		setwd(Jpath); command = paste("java -Xmx1600m MASTERCLASS  -f ", file, " -r ", finName, " -o ", finRep, " -abs ", absCut, sep=""); system(command);
		if(length(count.fields(finRep))>0) {
			re = read.table(finRep, sep="\t")
			for(z in 1:length(re[,1])) { re[z,4] = mean(r[re[z,6]:re[z,7],4]); }
			repdat= rbind(repdat, re)
		}
	}
	d = d[complete.cases(d),]; d[,1] = as.integer(d[,1]); d[,2] = as.integer(d[,2]); d[,3] = as.integer(d[,3]); d[,4] = as.double(d[,4]);
	write.table(data.frame(as.integer(d[,1]),as.integer(d[,2]),as.integer(d[,3]), d[,-(1:3)]), file=file, sep="\t", row.names=F, col.names=F, quote=F)	
	rep = rep[complete.cases(rep),]; rep[,1] = as.integer(rep[,1]); rep[,2] = as.integer(rep[,2]); rep[,3] = as.integer(rep[,3]); rep[,4] = as.double(rep[,4]);
	write.table(data.frame(as.integer(rep[,1]),as.integer(rep[,2]),as.integer(rep[,3]), rep[,-(1:3)]), file=finName, sep="\t", row.names=F, col.names=F, quote=F);
		if(length(repdat)>0) {
			repdat=repdat[complete.cases(repdat),]; repdat[,1] = as.integer(repdat[,1]); repdat[,2] = as.integer(repdat[,2]); repdat[,3] = as.integer(repdat[,3]); repdat[,4] = as.double(repdat[,4]);
			write.table(data.frame(as.integer(repdat[,1]),as.integer(repdat[,2]),as.integer(repdat[,3]), repdat[,-(1:3)]), file=finRep, sep="\t", row.names=F, col.names=F, quote=F)
		} else {
			write.table(NULL, file=finRep, sep="\t", row.names=F, col.names=F, quote=F)	
		}
	} else {
		write.table(NULL, file=finRep, sep="\t", row.names=F, col.names=F, quote=F)	
	
	}
}



#`makeFObs` <- function(file, idir, odir, conNum=2, methods = c("ADM2", "ADM3", "safe", "walk", "SMUG", "GADA", "SMAP", "dnacopy", "fsegment", "fast", "bcp", "CNCP"), absCut=0.3, gold = "../extdata/42Mcalls_all_feb09.txt") {
#
#	command = "cat";
#	rcommand = "rm";
#	
#	n = unlist(strsplit(file, "/"))
#	finName = paste(odir, substr(n[length(n)], 0, nchar(n[length(n)])-4), "_Allfeatures.txt", sep="")
#	finRep = paste(odir, substr(n[length(n)], 0, nchar(n[length(n)])-4), "_FinalReport.txt", sep="")
#
#	for(x in 1:length(methods)) {
#		
#		name1 = paste(odir, "N", substr(n[length(n)], 0, nchar(n[length(n)])-4), "_FinalReport_", methods[x], ".txt", sep="")
#		name2 = paste(odir, "P", substr(n[length(n)], 0, nchar(n[length(n)])-4), "_FinalReport_", methods[x], ".txt", sep="")
#		if(length(count.fields(name1))>0) {
#			d1 = read.table(name1, sep="\t")
#			d1=cbind(as.integer(d1[,1]), as.integer(d1[,2]), as.integer(d1[,3]), as.double(d1[,4]), as.character(methods[x]))
#			write.table(d1, file=name1, sep="\t", row.names=F, col.names=F, quote=F)
#		}
#		if(length(count.fields(name2))>0) {
#			d2 = read.table(name2, sep="\t")
#			d2=cbind(as.integer(d2[,1]), as.integer(d2[,2]), as.integer(d2[,3]), as.double(d2[,4]), as.character(methods[x]))
#			write.table(d2, file=name2, sep="\t", row.names=F, col.names=F, quote=F)
#		}
#		command = paste(command, name1, name2, sep=" ")
#		rcommand = paste(rcommand, name1, name2, sep=" ")
#	}
#	
#	command = paste(command, ">", finName, sep=" ")
#	system(command)
#	#system(rcommand)
#	
#	d = read.table(finName, sep="\t")
#	d[,1] = as.numeric(d[,1])
#	d=d[order(d[,1], d[,2], d[,3]),]
#	d=d[abs(d[,4])>absCut,]
#	d=d[complete.cases(d),]
#	write.table(d, file=finName, sep="\t", row.names=F, col.names=F, quote=F)
#	
#	n=unlist(strsplit(file, "/"))
#	n1 = paste(odir, "N", n[length(n)], sep="")
#	n2 = paste(odir, "P", n[length(n)], sep="")
#	tfile = paste(file, sep="")
#	command = paste("cat", n1, n2, ">", tfile, sep=" ")
#	system(command)
#	
#	#command = paste("rm", n1, n2, sep=" ")
#	#system(command)
#	
#	rep = read.table(finName, sep="\t")
#	d = read.table(tfile, sep="\t")
#	d=d[order(d[,1], d[,2], d[,3]),]
#	repdat = NULL
#	u = unique(rep[,1]) 
#	for(x in 1:length(u)) {	
#		r=d[d[,1]==u[x],]
#		rr = rep[rep[,1]==u[x],]
#		write.table(r, file=tfile, sep="\t", row.names=F, col.names=F, quote=F)
#		write.table(rr, file=finName, sep="\t", row.names=F, col.names=F, quote=F)
#		Jpath =  system.file("java", package="CNsolidate")
#		setwd(Jpath)
#		command = paste("java -Xmx1600m MASTERCLASS  -f ", tfile, " -r ", finName, " -o ", finRep, " -abs ", absCut, sep="")
#		system(command)
#		if(length(count.fields(finRep))>0) {
#			re = read.table(finRep, sep="\t")
#			for(z in 1:length(re[,1])) {
#				re[z,4] = mean(r[re[z,6]:re[z,7],4])
#			}
#			repdat= rbind(repdat, re)
#		}
#	}
#	d = d[complete.cases(d),]; d[,1] = as.integer(d[,1]); d[,2] = as.integer(d[,2]); d[,3] = as.integer(d[,3]); d[,4] = as.double(d[,4]);
#	rep = rep[complete.cases(rep),]; rep[,1] = as.integer(rep[,1]); rep[,2] = as.integer(rep[,2]); rep[,3] = as.integer(rep[,3]); rep[,4] = as.double(rep[,4]);
#	repdat=repdat[complete.cases(repdat),]; repdat[,1] = as.integer(repdat[,1]); repdat[,2] = as.integer(repdat[,2]); repdat[,3] = as.integer(repdat[,3]); repdat[,4] = as.double(repdat[,4]);
#	write.table(data.frame(as.integer(d[,1]),as.integer(d[,2]),as.integer(d[,3]), d[,-(1:3)]), file=tfile, sep="\t", row.names=F, col.names=F, quote=F)
#	write.table(data.frame(as.integer(rep[,1]),as.integer(rep[,2]),as.integer(rep[,3]), rep[,-(1:3)]), file=finName, sep="\t", row.names=F, col.names=F, quote=F)
#	write.table(data.frame(as.integer(repdat[,1]),as.integer(repdat[,2]),as.integer(repdat[,3]), repdat[,-(1:3)]), file=finRep, sep="\t", row.names=F, col.names=F, quote=F)
#
#}
