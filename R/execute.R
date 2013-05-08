`execute` <- function() UseMethod("execute")

`execute.default` <- function(file, set=NULL) {
	if(is.null(set)) {
		set = settings()
		set$files$file = file
	}
}

execute.normalise <- function(set=NULL) {
	if(is.null(set)) {
		set = settings()
	}
	if(set$files$file=="") {
		warning("No input file specificed!")
	}
	if(set$files$idir=="") {
		warning("No input directory specificed!")
	}
	if(set$files$odir=="") {
		warning("No output directory specificed!")
	}
	
	#do.call(names(set$norm[1]), list(set))
	
	for(x in 1:length(set$norm)) {
		if(set$norm[x] == 1) {
			do.call(names(set$norm[x]), list(set))
		}
	}
}

execute.algorithms <- function(set=NULL) {
	
	if(is.null(set)) {
		set = settings()
	}
	if(set$files$file=="") {
		warning("No input file specificed!")
	}
	if(set$files$idir=="") {
		warning("No input directory specificed!")
	}
	if(set$files$odir=="") {
		warning("No output directory specificed!")
	}
	
	for(x in 1:length(set$algorithms)) {
		if(set$algorithms[x] == 1) {
			result <- try( do.call(names(set$algorithms[x]), list(set)) )
			if(class(result) == "try-error")	{
				write.table(NULL, file=paste(substr(set$files$file, 1, nchar(set$files$file)-4), "_FinalReport_", names(set$algorithms[x]), ".txt", sep=""), row.names=F, col.names=F, quote=F)
				setwd(set$files$idir)
			}
		}
	}
}

execute.combine <- function(set=NULL) {
	
	if(is.null(set)) {
		set = settings()
	}
	if(set$files$file=="") {
		warning("No input file specificed!")
	}
	if(set$files$idir=="") {
		warning("No input directory specificed!")
	}
	if(set$files$odir=="") {
		warning("No output directory specificed!")
	}
	
	for(x in 1:length(set$combine)) {
		if(set$combine[x] == 1) {
			result <- try( do.call(names(set$combine[x]), list(set)) )
			if(class(result) == "try-error")	{
				write.table(NULL, file=paste(substr(set$files$file, 1, nchar(set$files$file)-4), "_FinalReport.txt", sep=""), row.names=F, col.names=F, quote=F)
				setwd(set$files$idir)
			}
		}
	}
	
}

         
execute.dddpipe <- function(set=NULL) {

	library(DataDB)
	set = NULL
	
	if(is.null(set)) {
		set = settings()
	}
	if(!class(set)=="settings") {
		 warning("setting object not of class settings!!!")
	}
	
	sanger_id = 1;
	
	# Get sample level info for proband - sanger_id
	aCGH_file1 = "/Users/tf2/Desktop/Packages/UKBS_nonNovel_data/pair/US09463735_253122010388_S01_CGH_105_Dec08.txt";
	aCGH_file2 = "/Users/tf2/Desktop/Packages/UKBS_nonNovel_data/pair/US09463735_253122110123_S01_CGH_105_Dec08.txt";
	
	 
	set$files$idir ="/Users/tf2/Desktop/Packages/UKBS_nonNovel_data/pair/"
	odir ="/Users/tf2/Desktop/Packages/UKBS_nonNovel_data/output/"
	
	outdir1 = paste(odir, "/2531220/", sep="")
	command = paste("mkdir ", outdir1, sep=""); system(command);
	command = paste("cp ", aCGH_file1, " ", outdir1, sep=""); system(command);
	n = unlist(strsplit(aCGH_file1, "/"))
	aCGH_file1 = paste(outdir1, "/", n[length(n)], sep="");
	
	# aCGH Dye Normalise, sample tracking and initial QC
	
	set$files$file = aCGH_file1; set$files$odir = outdir1;
	acghQC1(set$files$file,sanger_id)
	set$files$format = 'fe'; jspline(set); set$files$format='bed';
	data = inData(set$files$file, set$files$format)
	insert_chr_medians(set$files$file, sanger_id, chr_m(set$files$file, set$files$odir));
	acghQC2(set$files$file, sanger_id)
	split_array(set$files$file, set$files$odir)
	
	if(is.null(aCGH_file2)) {
		
		### Tracking single file 
		cnv_tracking_single(set$files$file, sanger_id)	
		
		} else {		
	
		outdir2 = paste(odir, "/2531221/", sep="")
		command = paste("mkdir ", outdir2, sep=""); system(command);
		command = paste("cp ", aCGH_file2, " ", outdir2, sep=""); system(command);
		n = unlist(strsplit(aCGH_file2, "/"))
		aCGH_file2 = paste(outdir2, "/", n[length(n)], sep="");
		
		set$files$file = aCGH_file2; set$files$odir = outdir2;
		acghQC1(set$files$file,sanger_id)
		set$files$format='fe'; jspline(set); set$files$format='bed';
		data=inData(set$files$file, set$files$format)
		insert_chr_medians(set$files$file, sanger_id, chr_m(set$files$file, set$files$odir));
		acghQC2(set$files$file, sanger_id)
		split_array(set$files$file, set$files$odir)
		
		### Tracking both files 
		cnv_tracking(aCGH_file1, aCGH_file2, sanger_id)	
	}
	
	# Self normalise
	
	# Slide1
	
	set$files$file = aCGH_file1; set$files$odir = outdir1;
	if( length(grep("2531220", aCGH_file1))>0 )  {
		set$norm$sfile = "";
	} else {
		set$norm$sfile = "";
	}
	#self(set);
	
	# Slide2
	set$files$file = aCGH_file2; set$files$odir = outdir2;
	if( length(grep("2531220", aCGH_file2))>0 )  {
		set$norm$sfile = "";
	} else {
		set$norm$sfile = "";
	}
	#self(set);
	
	
	# Wavelet Normalise
	set$files$file = aCGH_file1; set$files$odir = outdir1; wave(set);
	set$files$file = aCGH_file2; set$files$odir = outdir2; wave(set);
	
	# aCGH Discovery, Combination and scoring
	# Slide1
	set$files$file = aCGH_file1; set$files$odir = outdir1;
	execute.algorithms(set); execute.combine(set);
	# Slide2
	set$files$file = aCGH_file2; set$files$odir = outdir2;
	execute.algorithms(set); execute.combine(set);
	
	# De-novo estimatation
	# ......
	
	
	# Insert to Results DataDB
	#n = unlist(strsplit(aCGH_file1, "/"))
	#report1 = paste(outdir1, substr(n[length(n)], 0, nchar(n[length(n)])-4), "_FinalReport.txt", sep="")
	insertDetectionFile(aCGH_file1, sanger_id)
	
	#n = unlist(strsplit(aCGH_file2, "/"))
	#report2 = paste(outdir1, substr(n[length(n)], 0, nchar(n[length(n)])-4), "_FinalReport.txt", sep="")
	insertDetectionFile(aCGH_file2, sanger_id);
	
	# Update LIMS 
	# ...
	
	# Update Decipher
	# ...
}
