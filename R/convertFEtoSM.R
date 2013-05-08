`convertFEtoSM` <- function(dir) {
	d = getwd()
	Jpath =  system.file("java", package="CNsolidate")
	setwd(Jpath)
	command = paste("java -Xmx1600m FileConverter -d ", dir, sep="")
	system(command)
	setwd(d)
}