`exampleADM2` <- function(outfile="example_ouput_ADM2.txt") {
	exFile= paste(system.file("data", package = "CNsolidate"), "/", "example.bed", sep="");
	plotADM2(ADM2(file=exFile, outfile=outfile));
}