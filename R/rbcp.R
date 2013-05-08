`rbcp` <- 
function(file, idir, odir, reso=1000, pc=0.05, p1=0.95, p2=0.01, mp=5, mr=0.2, pt=3, gold="../data/42Mcalls_all_feb09.txt") {		
		
		ifile = paste(idir, file, sep="")
		m = rSegmentation(ifile, reso, 1)
		
		A = mergeAll(m, pc, p1, p2)
		t = A[[1]]
		for(a in 2:length(A)) {
			t = rbind(t, A[[a]])	
		}
		
		cur = getwd()
		name = paste(odir, substr(file, 0, nchar(file)-4), ".temp", sep="")
		t[,1] = as.integer(t[,1])
		t[,2] = as.integer(t[,2])
		t[,3] = as.integer(t[,3])
		write.table(t, file=name, sep="\t", row.names=F, col.names=F, quote=F)

		Jpath =  system.file("java", package="CNsolidate")
		setwd(Jpath)
		command = paste("/software/bin/java -Xmx1600m Clean -f ", name, " -mp ", mp, " -mr ", mr, " -pc ", pc, " -pt ", pt, " -gold ", gold, " -ty BCP", sep="")
		system(command)
		setwd(cur)
	}


