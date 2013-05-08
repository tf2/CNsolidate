`jspline` <- function(set=NULL, outfile=NULL) {
	t = paste(set$files$file, ".temp", sep=""); 
	cur = getwd(); 
	Jpath = system.file("java", package = "CNsolidate"); 
	setwd(Jpath)
	command =paste("java -Xmx1600m Spline -f ", set$files$file, " -t ", t, " -fo ", set$files$format, " -th ", set$norm.settings$sthes, " -fa ", set$norm.settings$sFact, " -kn ", set$norm.settings$sknot, " -it ", set$norm.settings$sIt, sep="" ); 
	system(command); command = paste("mv ", t, " ", outfile, sep=""); system(command); setwd(cur);
}

#`jspline` <- function(file, idir, format="bed", thres=0.68, fact=4.5, knots=1000, iter=5) {
#	t = paste(file, ".temp", sep=""); cur = getwd(); Jpath = system.file("java", package = "CNsolidate"); setwd(Jpath)
#	command =paste("java -Xmx1600m Spline -f ",file, " -t ", t, " -fo ", format, " -th ", thres, " -fa ", fact, " -kn ", knots, " -it ", iter, sep="" ); 
#	system(command); command = paste("mv ", t, file); system(command); setwd(cur);
#}
