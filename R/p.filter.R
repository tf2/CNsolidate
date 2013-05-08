`p.filter` <- 
function(idir, odir) {		
		
		cur = getwd()
		Ppath =  system.file("perl", package="CNsolidate")
		setwd(Ppath)
		command = paste("perl makeFilter.pl ", idir, " ", odir, sep="")
		system(command)
		setwd(cur)
	}


