`self` <- function(set=NULL) {
	data = set$data
	if(is.null(data)) { data = read.table(set$files$file); }
	if(!is.null(set$norm.setting$sfile)) {
		self=read.table(set$norm.setting$sfile);
		log2_dye = data[,4]
		data[,4]=residuals(lm(data[,4]~self[,1]))
		u =unique(data[,1]); dat = NULL;
		for(x in 1:length(u)) {
			me = median(data[data[,1]==u[x],4])
			# exclude PAR regions - hg19
			if(u[x]==23) {
				chrX = data[data[,1]==u[x],]
				chrX=chrX[chrX[,2]<60001 | chrX[,2]>2699520,]
				chrX=chrX[chrX[,2]<154931044 | chrX[,2]>155260560,4]
				me = median(chrX)
			}
			if(u[x]==24) {
				chrY = data[data[,1]==u[x],]
				chrY=chrY[chrY[,2]<10001 | chrY[,2]>2649520,]
				chrY=chrY[chrY[,2]<59034050 | chrY[,2]>59363566,4]
				me = median(chrY)
			}
			
			data[data[,1]==u[x],4] = data[data[,1]==u[x],4]-me
			#if(length(data[data[,1]==u[x],1])>5000) {
				dat=rbind(dat, data[data[,1]==u[x],])
			#}	
		}
		dat= cbind(dat[,1:5], log2_dye, dat[,-(1:5)])
		write.table(dat, file=set$files$file, sep="\t", row.names=F, col.names=F, quote=F)
	}	
}

#`self` <- function(data, file, self) {		
#		self=read.table(self);
#		data[,4]=residuals(lm(data[,4]~self[,1]))
#		u =unique(data[,1]); dat = NULL;
#		for(x in 1:length(u)) {
#			data[data[,1]==u[x],4] = data[data[,1]==u[x],4]-median(data[data[,1]==u[x],4])
#			if(length(data[data[,1]==u[x],1])>5000) {
#				dat=rbind(dat, data[data[,1]==u[x],])
#			}	
#		}
#		write.table(dat, file=file, sep="\t", row.names=F, col.names=F, quote=F)
#}
