`weight.algorithms2` <- function(file, odir) {	
	n = unlist(strsplit(file, "/"))
	report = paste(odir, substr(n[length(n)], 0, nchar(n[length(n)])-4), "_FinalReport.txt", sep="")
	numberMatrix = paste(odir, "/", "Cmatrices/", substr(n[length(n)], 0, nchar(n[length(n)])-4), "_NumberMatrix.txt", sep="")
	log = read.table(file); data = read.table(numberMatrix, header=T); data = data[order(colnames(data))]
	for(x in 1:length(data[1,])) {
		data[data[,x]>1,x] = 1
	}
	data = as.matrix(data); rp = quantile(abs(log[,4]), probs=0.68); dl = dLRs(log[,4]); wa = ws(log[,4])	
	cenVt = read.table("../data/central_NoiseWeights_rp68.txt", header=T)
	difVt = read.table("../data/difference_NoiseWeights_dLRs.txt", header=T)
	wavVt = read.table("../data/wave_NoiseWeights_ws.txt", header=T)
	cen = vector(); dif = vector(); wav = vector();
	for(y in 1:length(data[1,])) {
		cen[y]= cenVt[which.min(abs(rp-cenVt[,1])), colnames(cenVt) == colnames(data)[y]]
		dif[y]= difVt[which.min(abs(dl-difVt[,1])), colnames(difVt) == colnames(data)[y]]
		wav[y]= wavVt[which.min(abs(wa-wavVt[,1])), colnames(wavVt) == colnames(data)[y]]
	}
	
	score = ((data%*%cen) + (data%*%dif) + (data%*%wav))/3;
	score = score/max(score); r=cbind(read.table(report),score);
	write.table(r, file=report, sep="\t", row.names=F, col.names=F, quote=F)
}