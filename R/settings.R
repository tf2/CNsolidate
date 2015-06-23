`settings` <- function() UseMethod("settings")

`settings.default` <- function() {
	settings <- list (
	
	"files" =list("file" = "", "idir" = "", "odir"="", "format"= "bed", "gold"="../extdata/all_sets_merged_hg19.txt", "gapfile"="../extdata/centromeres_hg19.txt"),	
	
	"norm" = list("jspline" = 1, "self" = 1, "gc" = 0, "wave" = 1),
		"norm.settings" =list("sFact" = 4.5, "sthes" = 0.68, "sknot" = 1000, "sIt" = 5, "gfile" = NULL, "sfile" = NULL, "wFac" = 1.2, "wFixed" = 0),		
	
	"algorithms" = list("rbcp" = 0, "safe" = 0, "fast" = 0, "ADM3" = 1, "SMA" = 1, "GADA" = 1, "fsegment" = 0, "rcncp" = 0, "dnacopy" = 1, "SMUG" = 1, "SMAP" = 1, "walk" = 0),
		"algorithms.sep" = list("rbcp" = 0, "safe" = 1, "fast" = 1, "ADM3" = 0, "SMA" = 0, "GADA" = 0, "fsegment" = 1, "rcncp" = 1, "dnacopy" = 1, "SMUG" = 1, "SMAP" = 1, "walk" = 0),
		#"rbcp" = list(),  
		"safe" = list("mL" = 3, "mR" = 0.35),   				
		"fast" = list("mL" = 3, "mR" = 0.35, "pT" = 1), 
		"ADM3" = list("t" = 5, "np" = 3, "mR" = 0.35, "nd" = 3),
		"SMA" = list("s" = 5, "np" = 3, "mR" = 0.35, "fact" = 3),
		"GADA" = list("alpha" = 0.1, "T" = 3, "mL" = 3, "mR" = 0.35),
		"fsegment" = list("sL" = 1000, "cpP" = 2, "mL" = 3, "mR" = 0.35), 
		"rcncp" = list("fac" = 3, "a" = 3, "n" = 2, "mR" = 0.35, "prob" = 0.4, "dist" = 5000),  				
		"dnacopy" = list("smR" = 2, "oSD" = 4, "sSD" = 2, "trim" = 0.025, "mL" = 3, "mR" = 0.35),
		"SMUG" = list("reso" = 1000, "pc" = 0.05, "p1" = 0.95, "p2" = 0.01, "mL" = 3, "mR" = 0.35, "pT" = 3),  
		"SMAP" = list("reso" = 1000, "pc" = 0.05, "p1" = 0.95, "p2" = 0.01, "mL" = 3, "mR" = 0.35, "pT" = 3),  			   				  				
		"walk" = list("tr1" = 2, "len1"=25, "len2"=5, "wlen1"=100, "wlen2"=500, "siZ"=500, "mL" = 3, "pc" = 0.05, "mR" = 0.35, "pT" = 3  ),  	
   	
   	"combine" =list("makeFObs" = 1,  "mapBreakc" = 1, "gap" = 1, "overlap" = 1, "derive.weights" = 1, "localp" = 1, "combin.image" = 0, "makeGff" = 0, "QC" = 0, "overlap3" =0),
   		"combine.settings" =list("Consenues" = 2, "absRatio" = 0.3, "minProbe" = 3, "lRatio" = 0.4, "gRatio" = 0.5, "lest" = 250, "gest" = 200, "tfac" = 2, "adj" = 1, "level" = 400), 	
	
	"data" = NULL
	
	)
	class(settings) <- "settings";
invisible(settings)
}

`settings.v1` <- function() {
	settings <- list (
	
	"files" =list("file" = "", "idir" = "", "odir"="", "format"= "bed", "gold"="../extdata/all_sets_merged_hg19.txt", "gapfile"="../extdata/centromeres_hg19.txt"),	
	
	"norm" = list("jspline" = 1, "self" = 1, "gc" = 0, "wave" = 1),
		"norm.settings" =list("sFact" = 4.5, "sthes" = 0.68, "sknot" = 1000, "sIt" = 5, "gfile" = NULL, "sfile" = NULL, "wFac" = 1.2, "wFixed" = 0),		
	
	"algorithms" = list("rbcp" = 0, "safe" = 0, "fast" = 0, "ADM3" = 1, "SMA" = 1, "GADA" = 1, "fsegment" = 0, "rcncp" = 0, "dnacopy" = 1, "SMUG" = 1, "SMAP" = 1, "walk" = 0),
		"algorithms.sep" = list("rbcp" = 0, "safe" = 1, "fast" = 1, "ADM3" = 0, "SMA" = 0, "GADA" = 0, "fsegment" = 1, "rcncp" = 1, "dnacopy" = 1, "SMUG" = 1, "SMAP" = 1, "walk" = 0),
		#"rbcp" = list(),  
		"safe" = list("mL" = 3, "mR" = 0.35),   				
		"fast" = list("mL" = 3, "mR" = 0.35, "pT" = 1), 
		"ADM3" = list("t" = 5, "np" = 3, "mR" = 0.35, "nd" = 3),
		"SMA" = list("s" = 5, "np" = 3, "mR" = 0.35, "fact" = 3),
		"GADA" = list("alpha" = 0.1, "T" = 3, "mL" = 3, "mR" = 0.35),
		"fsegment" = list("sL" = 1000, "cpP" = 2, "mL" = 3, "mR" = 0.35), 
		"rcncp" = list("fac" = 3, "a" = 3, "n" = 2, "mR" = 0.35, "prob" = 0.4, "dist" = 5000),  				
		"dnacopy" = list("smR" = 2, "oSD" = 4, "sSD" = 2, "trim" = 0.025, "mL" = 3, "mR" = 0.35),
		"SMUG" = list("reso" = 1000, "pc" = 0.05, "p1" = 0.95, "p2" = 0.01, "mL" = 3, "mR" = 0.35, "pT" = 3),  
		"SMAP" = list("reso" = 1000, "pc" = 0.05, "p1" = 0.95, "p2" = 0.01, "mL" = 3, "mR" = 0.35, "pT" = 3),  			   				  				
		"walk" = list("tr1" = 2, "len1"=25, "len2"=5, "wlen1"=100, "wlen2"=500, "siZ"=500, "mL" = 3, "pc" = 0.05, "mR" = 0.35, "pT" = 3  ),  	
   	
#	"combine" =list("makeFObs" = 1,  "mapBreakc" = 1, "gap" = 1, "overlap" = 1, "weight.algorithms_across_samples" = 1, "localp" = 1, "combin.image" = 1, "makeGff" = 1, "QC" = 1, "overlap3" =1),

#	"combine" =list("makeFObs" = 1,  "mapBreakc" = 1, "gap" = 1, "overlap" = 1, "weight.algorithms_across_samples" = 1, "localp" = 1, "CI.breakpoints" =1, "freq.overlap" =1),
					  
	"combine" =list("makeFObs" = 1,  "mapBreakc" = 1, "gap.depends" = 1, "overlap" = 1, "weight.algorithms_across_samples" = 1, "localp" = 1, "CI.breakpoints" =1, "freq.overlap" =1),
					  
		"combine.settings" =list("Consenues" = 2, "absRatio" = 0.3, "minProbe" = 3, "lRatio" = 0.4, "gRatio" = 0.5, "lest" = 250, "gest" = 200, "tfac" = 2, "adj" = 1, "level" = 400), 	
	
	"data" = NULL
	
	)
	class(settings) <- "settings";
invisible(settings)
}


`settings.independant` <- function() {
	settings <- list (
					  
					  "files" =list("file" = "", "idir" = "", "odir"="", "format"= "bed", "gold"="../extdata/all_sets_merged_hg19.txt", "gapfile"="../extdata/centromeres_hg19.txt"),	
					  
					  "norm" = list("jspline" = 1, "self" = 1, "gc" = 0, "wave" = 1),
					  "norm.settings" =list("sFact" = 4.5, "sthes" = 0.68, "sknot" = 1000, "sIt" = 5, "gfile" = NULL, "sfile" = NULL, "wFac" = 1.2, "wFixed" = 0),		
					  
					  "algorithms" = list("rbcp" = 0, "safe" = 0, "fast" = 0, "ADM3" = 1, "SMA" = 1, "GADA" = 1, "fsegment" = 0, "rcncp" = 0, "dnacopy" = 1, "SMUG" = 1, "SMAP" = 1, "walk" = 0),
					  "algorithms.sep" = list("rbcp" = 0, "safe" = 1, "fast" = 1, "ADM3" = 0, "SMA" = 0, "GADA" = 0, "fsegment" = 1, "rcncp" = 1, "dnacopy" = 1, "SMUG" = 1, "SMAP" = 1, "walk" = 0),
#"rbcp" = list(),  
					  "safe" = list("mL" = 3, "mR" = 0.35),   				
					  "fast" = list("mL" = 3, "mR" = 0.35, "pT" = 1), 
					  "ADM3" = list("t" = 5, "np" = 3, "mR" = 0.35, "nd" = 3),
					  "SMA" = list("s" = 5, "np" = 3, "mR" = 0.35, "fact" = 3),
					  "GADA" = list("alpha" = 0.1, "T" = 3, "mL" = 3, "mR" = 0.35),
					  "fsegment" = list("sL" = 1000, "cpP" = 2, "mL" = 3, "mR" = 0.35), 
					  "rcncp" = list("fac" = 3, "a" = 3, "n" = 2, "mR" = 0.35, "prob" = 0.4, "dist" = 5000),  				
					  "dnacopy" = list("smR" = 2, "oSD" = 4, "sSD" = 2, "trim" = 0.025, "mL" = 3, "mR" = 0.35),
					  "SMUG" = list("reso" = 1000, "pc" = 0.05, "p1" = 0.95, "p2" = 0.01, "mL" = 3, "mR" = 0.35, "pT" = 3),  
					  "SMAP" = list("reso" = 1000, "pc" = 0.05, "p1" = 0.95, "p2" = 0.01, "mL" = 3, "mR" = 0.35, "pT" = 3),  			   				  				
					  "walk" = list("tr1" = 2, "len1"=25, "len2"=5, "wlen1"=100, "wlen2"=500, "siZ"=500, "mL" = 3, "pc" = 0.05, "mR" = 0.35, "pT" = 3  ),  	
					  
					  "combine" =list("makeFObs" = 1,  "mapBreakc" = 1, "gap" = 1, "overlap" = 1, "derive.weights.independant" = 1, "localp" = 1, "combin.image" = 1, "makeGff" = 1, "QC" = 1, "overlap3" =1),
					  "combine.settings" =list("Consenues" = 2, "absRatio" = 0.3, "minProbe" = 3, "lRatio" = 0.4, "gRatio" = 0.5, "lest" = 250, "gest" = 200, "tfac" = 2, "adj" = 1, "level" = 400), 	
					  
					  "data" = NULL
					  
					  )
	class(settings) <- "settings";
	invisible(settings)
}
