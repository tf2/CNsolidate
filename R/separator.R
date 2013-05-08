`separator` <-
function(file, idir, odir, ConNum=3, gold="../data/42Mcalls_all_feb09.txt") {		
  			 
	cur = getwd()
	#inFiles = dir(idir)
	#inFiles = inFiles[grep("N", inFiles)]
	
	#for( x in 1:length(inFiles)) {
	
		rep = paste(substr(file, 0, nchar(file)-4), "_FinalReport_Consensus_Report.txt", sep="") 
		inF = paste(idir, "/", file, sep="")
		outF = paste(odir, "/", rep, sep="")
		Jpath =  system.file("java", package="CNsolidate")
		setwd(Jpath)
		
		command = paste("java -Xmx1600m Separator -f ", inF, " -r ", outF, sep="")
		system(command)
		setwd(cur)
	#}
}
