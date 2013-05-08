`overlap2` <- 
function(file, idir, odir,methods = c("safe", "walk", "SMUG", "GADA", "SMAP", "dnacopy", "fsegment", "fast", "bcp", "CNCP"), gold = "../extdata/42Mcalls_all_feb09.txt") {		

		n = unlist(strsplit(file, "/"))
		file1 = paste(odir, substr(n[length(n)], 0, nchar(n[length(n)])-4), "_Allfeatures.txt", sep="")
		file2 = paste(odir, substr(n[length(n)], 0, nchar(n[length(n)])-4), "_FinalReport.txt", sep="")
		alf = paste(odir, substr(n[length(n)], 0, nchar(n[length(n)])-4), "_AllFDR.txt", sep="")
		rmc = "rm "
		catc = "cat "
		for(x in 1:length(methods)) {
		
			name = paste(odir, "/", substr(n[length(n)], 0, nchar(n[length(n)])-4), "_FinalReport_", methods[x], ".txt", sep="")
			command = paste("grep -i ", methods[x], " ",  file1, ">", name, sep="")
			system(command)
			Jpath =  system.file("java", package="CNsolidate")
			setwd(Jpath)
			command = paste("java -Xmx1600m Overlap2 -f ", name, " -o ",  name, " -al ", methods[x], " -g ", gold, sep="") 
			system(command)
			rmc = paste(rmc, name, " ", sep="")
			#system(command)
			catc = paste(catc, name, " ", sep="")
			#system(command)
			
		}
			command = paste("java -Xmx1600m Overlap2 -f ", file2, " -o ",  "combineFDR.txt", " -al ", "combined", " -g ", gold, sep="") 
			system(command)
			catc = paste(catc, "combineFDR.txt", " ", sep="")
			rmc = paste(rmc, "combineFDR.txt", " ", sep="")
		catc = paste(catc, ">", alf, sep="")
	system(catc)
	system(rmc)
}
