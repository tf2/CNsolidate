`split_probes` <- function(file, commonProbes = "../data/031220_031221_wtccc_and_sequenom_probe_ids.txt") {
	cur = getwd();
	Ppath =  system.file("perl", package="CNsolidate")
	setwd(Ppath)
		command = paste("perl split_common_rare_probe.pl ", file, " ", commonProbes, sep="")
		system(command)
	setwd(cur)
}