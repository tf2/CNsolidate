`dt.g.wavelet` <-
function(data, file, idir, odir, segMethod = "GADA") {		
  			 
	cur = getwd()
	name = NULL
	
	if(segMethod=="GADA") {
		GADA(data, file, idir, odir)
		n = unlist(strsplit(file, "/"))
		name = paste(odir, substr(n[length(n)], 0, nchar(n[length(n)])-4), "_FinalReport_GADA.txt", sep="")
	}
	
	Jpath =  system.file("java", package="CNsolidate")
	setwd(Jpath)
	
	command = paste("java -Xmx1600m genWave -f ", file, " -r ", name, sep="")

	system(command)
	setwd(cur)
}

