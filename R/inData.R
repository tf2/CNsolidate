`inData` <- function(file, format="gff", med_norm="global") {

	dat = NULL;

	if(format=="gff") {
	
		data = read.table(file)
		data=data[data[,1]!="chrXD",]
		data=data[data[,1]!="chr2L",]
		c = substr(data[,1], 4, 10)
		c[c=="X"] = 23
		c[c=="Y"] = 24
		c = as.numeric(c)
		data[,1] = c
		data = data[order(data[,1], data[,4], data[,5]),]
		data =cbind(data[,1], data[,4:6])
		data[,1] = as.integer(data[,1])
		data[,2] = as.integer(data[,2])
		data[,3] = as.integer(data[,3])
		data=data[complete.cases(data),]
		u =unique(data[,1]); dat = NULL;
		if(med_norm == "global") {
			dat=data
			dat[,4] = data[,4]-median(data[data[,1]<23,4])
		}
		if(med_norm == "chromosome") {
			for(x in 1:length(u)) {
				data[data[,1]==u[x],4] = data[data[,1]==u[x],4]-median(data[data[,1]==u[x],4])
				dat=rbind(dat, data[data[,1]==u[x],])
			}
		}
			write.table(dat, file=file, sep="\t", row.names=F, col.names=F, quote=F)
	}
	
	if(format=="bed" | format=="fe") {
		data = read.table(file)
		data=data[data[,1]!="chrXD",]
		data=data[data[,1]!="chr2L",]
		c = substr(data[,1], 4, 10)
		c[c=="X"] = 23
		c[c=="Y"] = 24
		c = as.numeric(c)
		data[,1] = c
		data = data[order(data[,1], data[,2], data[,3]),]
		data[,1] = as.integer(data[,1])
		data[,2] = as.integer(data[,2])
		data[,3] = as.integer(data[,3])
		data=data[complete.cases(data),]
		u =unique(data[,1]); dat = NULL;
		if(med_norm == "global") {
			dat=data
			dat[,4] = data[,4]-median(data[data[,1]<23,4])
		}
		if(med_norm == "chromosome") {
			for(x in 1:length(u)) {
				data[data[,1]==u[x],4] = data[data[,1]==u[x],4]-median(data[data[,1]==u[x],4])
				dat=rbind(dat, data[data[,1]==u[x],])
			}
		}
		write.table(dat, file=file, sep="\t", row.names=F, col.names=F, quote=F)
	}
	
	invisible(dat)
}