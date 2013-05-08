`CI.breakpoints` <- function(set=NULL, res=5) {
	n = unlist(strsplit(set$files$file, "/"));
	report = paste(set$files$odir, "/", substr(n[length(n)], 0, nchar(n[length(n)])-4), "_FinalReport.txt", sep="")
	dd<-read.table(set$files$file); rep<-read.table(report);
	d = dd[,1:4]; pin = 1; u = unique(rep[,1])
	mat = matrix(nrow=length(rep[,1]), ncol=2)
	for(k in 1:length(u)) {
		chr<-u[k]
		mybpt<-rep[rep[,1]==chr,]
		x=d[d[,1]==chr,]
		names(x)<-c("chr","physposStart", "physposEnd", "residual")

		for(i in 1:length(mybpt[,1])) {

		lstart<-mybpt[i,2]
		lbpts<-make.bpts(x, lstart, res, 2)

		llog.lik.bpts<-vector()
		for (j in 1:length(lbpts)){
  			ldat<-x[x$physposStart>lbpts[1] & x$physposStart<lbpts[j],4]
  			rdat<-x[x$physposStart<lbpts[length(lbpts)] & x$physposStart>lbpts[j],4]
  			llog.lik.bpts[j]<-break.log.liks2(ldat,rdat)
		}

		lbpts<-lbpts[is.na(llog.lik.bpts)==F]
		llog.lik.bpts<-llog.lik.bpts[is.na(llog.lik.bpts)==F]
		lpprob<-bayes.break(llog.lik.bpts)
	
		rstart<-mybpt[i,3]
		rbpts<-make.bpts(x,rstart,res, 3)

		log.lik.bpts<-vector()
		for (j in 1:length(rbpts)){
  			ldat<-x[x$physposEnd>rbpts[1] & x$physposEnd<rbpts[j],4]
  			rdat<-x[x$physposEnd<rbpts[length(rbpts)] & x$physposEnd>rbpts[j],4]
  			log.lik.bpts[j]<-break.log.liks2(ldat,rdat)
		}

		rbpts<-rbpts[is.na(log.lik.bpts)==F]
		log.lik.bpts<-log.lik.bpts[is.na(log.lik.bpts)==F]
		rpprob<-bayes.break(log.lik.bpts)
	
		my.map<-x[x$physposStart> lstart-20000 & x$physposEnd <= rstart+20000,]
	
		mat[pin,1]<-break.CI(llog.lik.bpts,lpprob,lbpts,my.map,1)
		mat[pin,2]<-break.CI(log.lik.bpts,rpprob,rbpts,my.map,2)
		pin=pin+1
		}
	}
	ndat = cbind(rep,mat);
	write.table(ndat, file=report, sep="\t", row.names=F, col.names=F, quote=F)
invisible(ndat)
}