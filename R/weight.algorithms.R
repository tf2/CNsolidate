`weight.algorithms` <- function(file, odir) {
	
	n = unlist(strsplit(file, "/"))
	report = paste(odir, substr(n[length(n)], 0, nchar(n[length(n)])-4), "_FinalReport.txt", sep="")
	numberMatrix = paste(odir, "/", "Cmatrices/", substr(n[length(n)], 0, nchar(n[length(n)])-4), "_NumberMatrix.txt", sep="")
	
	data = read.table(numberMatrix, header=T)
	data = data[order(colnames(data))]
	for(x in 1:length(data[1,])) {
		data[data[,x]>1,x] = 1
	}

	score = vector()
	wm = read.table(paste(system.file("data", package="CNsolidate"), "/", "WeightMatrix.txt", sep=""), row.names=1, header=T)
	vt = read.table(paste(system.file("data", package="CNsolidate"), "/", "vectorLookUpTable.txt", sep=""), header=T)
	
	rep = read.table(report)
	log = read.table(file)
	noi = dLRs(log[,4])

	nv = vector()
	pin=2
	e = length(vt[1,])-1
	index = which.min(abs(noi-vt[,1]))
		
		for(x in 1:e) {
			nv[x] = vt[index,pin]
			pin=pin+1
		}

	nwm = nv*wm
		for(pin in 1:length(data[,1])) {
			s=0
			count = sum(data[pin,])
			for(x in 1:length(nwm[1,])) {
				s = s+ sum(nwm[x,]*data[pin,])
			}

		#score[pin] = s
		score[pin] = s/count
		}
		
	#rep = cbind(rep, score/max(score))
	rep = cbind(rep, score)
	write.table(rep, file=report, sep="\t", row.names=F, col.names=F, quote=F)
	
return(rep)
}




