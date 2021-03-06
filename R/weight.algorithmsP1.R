`weight.algorithmsP1` <- function(file, odir) {	
	n = unlist(strsplit(file, "/"))
	report = paste(odir, substr(n[length(n)], 0, nchar(n[length(n)])-4), "_FinalReport.txt", sep="")
	numberMatrix = paste(odir, "/", "Cmatrices/", substr(n[length(n)], 0, nchar(n[length(n)])-4), "_NumberMatrix.txt", sep="")
	log = read.table(file); data = read.table(numberMatrix, header=T); data = data[order(colnames(data))]
	for(x in 1:length(data[1,])) {
		data[data[,x]>1,x] = 1
	}
	data = as.matrix(data); rp = quantile(abs(log[,4]), probs=0.68); dl = dLRs(log[,4]); wa = quantile(abs(runmed(log[,4], 301)), probs=0.68)	
	cen = vector(); dif = vector(); wav = vector();
	for(y in 1:length(data[1,])) {
		cen[y]= centralNoiseWeights(rp, colnames(data)[y])
		dif[y]= differenceNoiseWeights(dl, colnames(data)[y])
		wav[y]= WaveNoiseWeights(wa, colnames(data)[y])
	}
	als = vector();
	for(x in 1:length(data[,1])) {
		als[x]=sum(data[1,]);
	} 
	als[als==0]=1;
	#score = ((data%*%cen) + (data%*%dif) + (data%*%wav) + (als)/length(data[1,]))/4;
	score = ( (data%*%cen) + (data%*%dif) + (data%*%wav) )/3;
	rep = read.table(report, sep="\t"); nsc = vector();
	for(x in 1:length(rep[,1])) {
		wi=getSizeWeight(rep[x,5])
		ra=getRatWeight(abs(rep[x,4]), dl)
		nsc[x]=(wi+ra)/2
	}
	score = (score+nsc)/2; 
	score = rangescale(score);
	r=cbind(rep,score);
	write.table(r, file=report, sep="\t", row.names=F, col.names=F, quote=F)
}