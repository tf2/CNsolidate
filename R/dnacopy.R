`dnacopy` <- function(set=NULL) {
  	data = set$data;
	if(is.null(data)) { data = read.table(set$files$file); }
  		
  		dnac <- function(data, set) {
  			CNA.object <- CNA(cbind(data[,4]), paste("chr", data[,1], sep=""), data[,2], data.type="logratio",sampleid="sam1" )	
  			smoothed.CNA.object <- smooth.CNA(CNA.object)
  			segment.smoothed.CNA.object <- segment(smoothed.CNA.object, verbose=1)
			rep = segment.smoothed.CNA.object$output; rep = cbind(rep[,2:4], rep[,6], rep[,5]);
			rep = rep[abs(rep[,4])>set$dnacopy$mR,]; rep = rep[rep[,5]>set$dnacopy$mL,]
			rep[,1] = substr(rep[,1], 4 , 10); rep[,1] = as.integer(rep[,1]);
			rep[,2] = as.integer(rep[,2]); rep[,3] = as.integer(rep[,3]);
		invisible(rep)
		}
	
	n = unlist(strsplit(set$files$file, "/"))
	name = paste(set$files$odir, "/", substr(n[length(n)], 0, nchar(n[length(n)])-4), "_FinalReport_dnacopy.txt", sep="")
	
	rep=NULL;
	if(set$algorithms.sep$dnacopy==1) {
		rep = rbind(dnac(data[data[,4]>0,], set), dnac(data[data[,4]<0,], set)); 
	} else {
		rep = dnac(data[data[,4]>0,], set); 
	}
	if(length(rep[,1])>0) {
		rep=rep[order(rep[,1],rep[,2],rep[,3]),];
	}
	write.table(rep, file=name, sep="\t", quote=F, col.names=F, row.names=F)
  	invisible(rep)
  }


#`dna.copy` <- function(data, file, idir, odir, minL = 5, minR=0.2) {  	
#  	name = substr(file, 0 , nchar(file)-4)
#  	CNA.object <- CNA(cbind(data[,4]), paste("chr", data[,1], sep=""), data[,2], data.type="logratio",sampleid="sam1" )	
#  	smoothed.CNA.object <- smooth.CNA(CNA.object)
#  	segment.smoothed.CNA.object <- segment(smoothed.CNA.object, verbose=1)  	
#	rep = segment.smoothed.CNA.object$output
#	rep = cbind(rep[,2:4], rep[,6], rep[,5])
#	rep = rep[abs(rep[,4])>minR,]
#	rep = rep[rep[,5]>minL,]
#	rep[,1] = substr(rep[,1], 4 , 10)
#	rep[,1] = as.integer(rep[,1])
#	rep[,2] = as.integer(rep[,2])
#	rep[,3] = as.integer(rep[,3])
#	write.table(rep, file=paste(name, "_FinalReport_dnacopy.txt", sep=""), sep="\t", quote=F, col.names=F, row.names=F)
# 	return(rep)
# }
  