`norm_split_dir` <- function(idir, format="fe", wFac=1.2) {
	odir = idir
	files=dir(idir)
	for(x in 1:length(files)) {
		file=paste(idir, "/", files[x], sep="")
		jspline(file, idir, format)
		data = inData(file, "bed")
		wave(file, idir, odir, wFac)
		split_probes(file)
	}
}