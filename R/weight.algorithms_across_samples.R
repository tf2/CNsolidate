`weight.algorithms_across_samples` <- function(set = NULL) {	
	n = unlist(strsplit(set$files$file, "/"))
	report = paste(set$files$odir, "/", substr(n[length(n)], 0, nchar(n[length(n)])-4), "_FinalReport.txt", sep="")
	if(length(count.fields(report))>0) {
	numberMatrix = paste(set$files$odir, "/", "Cmatrices/", substr(n[length(n)], 0, nchar(n[length(n)])-4), "_NumberMatrix.txt", sep="")
	log = read.table(set$files$file); data = read.table(numberMatrix, header=T); data = data[order(colnames(data))]
	for(x in 1:length(data[1,])) {
		data[data[,x]>1,x] = 1
	}
	data = as.matrix(data); rp = quantile(abs(log[,4]), probs=0.68); dl = dLRs(log[,4]); 
	sp = 301
	if(sp>length(log[,4])) { sp = length(log[,4]/5); }
	wa = quantile(abs(runmed(log[,4], sp)), probs=0.68)	
	cen = vector(); dif = vector(); wav = vector();
	for(y in 1:length(data[1,])) {
		cen[y]= centralNoiseWeights(rp, colnames(data)[y])
		dif[y]= differenceNoiseWeights(dl, colnames(data)[y])
		wav[y]= WaveNoiseWeights(wa, colnames(data)[y])
	}
	als = vector();
	for(x in 1:length(data[,1])) {
		als[x]=sum(data[x,]);
	} 
	als[als==0]=1;
	
	#score = ((data%*%cen) + (data%*%dif) + (data%*%wav) + (als)/length(data[1,]))/4;
	score = vector();
	msc = matrix(ncol=4, nrow=length(data[,1]))
	for(x in 1 :length(data[,1])) {
		su = sum(data[x,])
		#score[x] = ( (sum(cen)/su) + (sum(dif)/su) + (sum(wav)/su) + (als[x])/length(data[1,]) ) /4
		msc[x,1] = (sum(cen*data[x,])/su)
		msc[x,2] = (sum(dif*data[x,])/su)
		msc[x,3] = (sum(wav*data[x,])/su)
		msc[x,4] = (als[x])/length(data[1,])
		score[x]=(msc[x,1]+msc[x,2]+msc[x,3]+msc[x,4])/4
	} 
	score[is.na(score)]=min(score,na.rm=T);
	
	#score = ( (data%*%cen) + (data%*%dif) + (data%*%wav) )/3;
	rep = read.table(report, sep="\t"); nsc = vector();
	for(x in 1:length(rep[,1])) {
		wi=getSizeWeight(rep[x,5])
		ra=getRatWeight(abs(rep[x,4]))
		nsc[x]=(wi+ra)/2
	}
	score = (score+nsc)/2; #score = ((1-dl)+score)/2; 
	score[is.na(score)]=min(score,na.rm=T);
	r=cbind(rep,score);
	if(set$combine.settings$adj==TRUE) {
		p = adj.f.t.n(dl,set$combine.settings$level)
		r[,length(r[1,])] = score-p
		p = adj.f.t.n.v2(dl)
		r = cbind(r,score-p)
		#r = cbind(r, adj.score(score,dl));
	}
	write.table(r, file=report, sep="\t", row.names=F, col.names=F, quote=F)
	}
}



#`weight.algorithms_across_samples` <- function(file, odir, adj=TRUE, level=400) {	
#
#	n = unlist(strsplit(file, "/"))
#	report = paste(odir, substr(n[length(n)], 0, nchar(n[length(n)])-4), "_FinalReport.txt", sep="")
#	numberMatrix = paste(odir, "/", "Cmatrices/", substr(n[length(n)], 0, nchar(n[length(n)])-4), "_NumberMatrix.txt", sep="")
#	log = read.table(file); data = read.table(numberMatrix, header=T); data = data[order(colnames(data))]
#	for(x in 1:length(data[1,])) {
#		data[data[,x]>1,x] = 1
#	}
#	data = as.matrix(data); rp = quantile(abs(log[,4]), probs=0.68); dl = dLRs(log[,4]); wa = quantile(abs(runmed(log[,4], 301)), probs=0.68)	
#	cen = vector(); dif = vector(); wav = vector();
#	for(y in 1:length(data[1,])) {
#		cen[y]= centralNoiseWeights(rp, colnames(data)[y])
#		dif[y]= differenceNoiseWeights(dl, colnames(data)[y])
#		wav[y]= WaveNoiseWeights(wa, colnames(data)[y])
#	}
#	als = vector();
#	for(x in 1:length(data[,1])) {
#		als[x]=sum(data[x,]);
#	} 
#	als[als==0]=1;
#	
#	#score = ((data%*%cen) + (data%*%dif) + (data%*%wav) + (als)/length(data[1,]))/4;
#	score = vector();
#	msc = matrix(ncol=4, nrow=length(data[,1]))
#	for(x in 1 :length(data[,1])) {
#		su = sum(data[x,])
#		#score[x] = ( (sum(cen)/su) + (sum(dif)/su) + (sum(wav)/su) + (als[x])/length(data[1,]) ) /4
#		msc[x,1] = (sum(cen*data[x,])/su)
#		msc[x,2] = (sum(dif*data[x,])/su)
#		msc[x,3] = (sum(wav*data[x,])/su)
#		msc[x,4] = (als[x])/length(data[1,])
#		score[x]=(msc[x,1]+msc[x,2]+msc[x,3]+msc[x,4])/4
#	} 
#	score[is.na(score)]=min(score,na.rm=T);
#	
#	#score = ( (data%*%cen) + (data%*%dif) + (data%*%wav) )/3;
#	rep = read.table(report, sep="\t"); nsc = vector();
#	for(x in 1:length(rep[,1])) {
#		wi=getSizeWeight(rep[x,5])
#		ra=getRatWeight(abs(rep[x,4]), dl)
#		nsc[x]=(wi+ra)/2
#	}
#	score = (score+nsc)/2; #score = ((1-dl)+score)/2; 
#	score[is.na(score)]=min(score,na.rm=T);
#	r=cbind(rep,score);
#	if(adj==TRUE) {
#		p = adj.f.t.n(dl,level)
#		r = cbind(r,score-p)
#		#r = cbind(r, adj.score(score,dl));
#	}
#	write.table(r, file=report, sep="\t", row.names=F, col.names=F, quote=F)
#}