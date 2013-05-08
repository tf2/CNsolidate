`makeGff` <- 
function(set=NULL, probe_score = "/nfs/ddd0/Tom/CODE/probe_scores/") {		
		file=set$files$file; odir=set$files$odir;
		cur = getwd()
		Jpath =  system.file("perl", package="CNsolidate")
		setwd(Jpath)
		n = unlist(strsplit(file, "/"))
		name = paste(odir, "/", n[length(n)], sep="")
		command = paste("perl probe_scoring.pl ", name, " ", probe_score, sep="")
		system(command)
		report = paste(odir, "/", substr(n[length(n)], 0, nchar(n[length(n)])-4), "_FinalReport.txt", sep="")
		command = paste("perl make_scores_GFF.pl ", name, " ", report, sep="")
		system(command)
		setwd(cur)
		d = read.table(file=set$files$file); d=d[,-(length(d[1,]))]; write.table(d, file=set$files$file, sep="\t", row.names=F, col.names=F, quote=F);
	}


