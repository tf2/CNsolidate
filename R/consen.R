`consen` <-
function(file, idir, odir, ConNum=3, gold="../data/42Mcalls_all_feb09.txt") {		
  			 
	cur = getwd()
	fin = paste(idir, file, sep="")
	name1 = paste(odir, substr(file, 0, nchar(file)-4), "_FinalReport_safe.txt", sep="")
	name2 = paste(odir, substr(file, 0, nchar(file)-4), "_FinalReport_GADA.txt", sep="")
	name3 = paste(odir, substr(file, 0, nchar(file)-4), "_FinalReport_SMUG.txt", sep="")
	name4 = paste(odir, substr(file, 0, nchar(file)-4), "_FinalReport_SMAP.txt", sep="")
	name5 = paste(odir, substr(file, 0, nchar(file)-4), "_FinalReport_walk.txt", sep="")

	Jpath =  system.file("java", package="CNsolidate")
	setwd(Jpath)
	#command = paste("java -Xmx1600m Over -f1 ", name1, " -f2 ", name2, " -f3 ", name3, sep="") #" -gold ", gold, sep="")

	command = paste("java -Xmx1600m mIntersect -f1 ", name1, " -f2 ", name2, " -f3 ", name3, " -f4 ", name4, " -f5 ", name5, " -cn ", ConNum, " -fin ", fin, " -gold ", gold, sep="")
	system(command)
	setwd(cur)
}

