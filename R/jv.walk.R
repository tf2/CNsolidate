`jv.walk` <- 
function(file, idir, odir, reso=1000, pc=0.05, p1=0.95, p2=0.01, mp=5, mr=0.2, pt=3, gold="../data/42Mcalls_all_feb09.txt") {		
		
		cur = getwd()
		Jpath =  system.file("java", package="CNsolidate")
		setwd(Jpath)
		command = paste("java -Xmx1600m ControlT -f ", file, " -mr ", mr, " -i ", idir, " -o ", odir, sep="")
		system(command)
		n = unlist(strsplit(file, "/"))
		name = paste(odir, substr(n[length(n)], 0, nchar(n[length(n)])-4), "_walk.temp", sep="")
		command = paste("java -Xmx1600m Clean -f ", name, " -mp ", mp, " -mr ", mr, " -pc ", pc, " -pt ", pt, " -gold ", gold, " -ty walk", sep="")
		system(command)
		setwd(cur)
	}


