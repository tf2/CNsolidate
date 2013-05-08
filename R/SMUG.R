`SMUG` <- function(set=NULL) {		
  	data = set$data
	if(is.null(data)) { data = read.table(set$files$file); }
	n = unlist(strsplit(set$files$file, "/"))
	name = paste(set$files$odir, "/", substr(n[length(n)], 0, nchar(n[length(n)])-4), "_FinalReport_SMUG.txt", sep="")
    smug <- function(data, set, name) {
    	m = rSegmentation(data, set$SMUG$reso)		 
    	A = mergeAll(m, set$SMUG$pc, set$SMUG$p1, set$SMUG$p2)
    	t = A[[1]]
    		if(length(A) > 1) {
				for(a in 2:length(A)) {
					t = rbind(t, A[[a]])	
				}
			}	 
		cur = getwd();
		t[,1] = as.integer(t[,1]); t[,2] = as.integer(t[,2]); t[,3] = as.integer(t[,3]);
		write.table(t, file=name, sep="\t", row.names=F, col.names=F, quote=F)
	
		Jpath =  system.file("java", package="CNsolidate"); setwd(Jpath);
		command = paste("java -Xmx1600m Clean -f ", name, " -mp ", set$SMUG$mL, " -mr ", set$SMUG$mR, " -pc ", set$SMUG$pc, " -pt ", set$SMUG$pT, " -gold ", set$files$gold, " -ty SMUG", sep="")
		system(command); setwd(cur);
	}
	
	r = NULL;
	if(set$algorithms.sep$SMUG==1) {
		smug(data[data[,4]>0,],set, name); 
		r1 = NULL; r2 = NULL;
		if(length(count.fields(name))>0) {
			r1= read.table(name);
		}
		smug(data[data[,4]<0,],set, name); 
		if(length(count.fields(name))>0) {
			r2= read.table(name);
		}
			r = rbind(r1,r2); 
	} else {
		smug(data,set, name);
		if(length(count.fields(name))>0) {
			r= read.table(name);
		}
	}
	if(length(r[,1])>0) {
		r=r[order(r[,1],r[,2],r[,3]),];
	}
	write.table(r, file=name, sep="\t", row.names=F, col.names=F, quote=F)
}


#`SMUG` <- function(data, file, idir, odir, reso=1000, pc=0.2, p1=0.95, p2=0.01, mp=5, mr=0.2, pt=3, gold="../data/42Mcalls_all_feb09.txt") {		
#  	
#    m = rSegmentation(data, reso)
#		 
#    A = mergeAll(m, pc, p1, p2)
#    t = A[[1]]
#    if(length(A) > 1) {
#		for(a in 2:length(A)) {
#			t = rbind(t, A[[a]])	
#		}
#	}	 
#	cur = getwd()
#	n = unlist(strsplit(file, "/"))
#	name = paste(odir, substr(n[length(n)], 0, nchar(n[length(n)])-4), "_SMUG.temp", sep="")
#	t[,1] = as.integer(t[,1])
#	t[,2] = as.integer(t[,2])
#	t[,3] = as.integer(t[,3])
#	write.table(t, file=name, sep="\t", row.names=F, col.names=F, quote=F)
#
#	Jpath =  system.file("java", package="CNsolidate")
#	setwd(Jpath)
#	command = paste("java -Xmx1600m Clean -f ", name, " -mp ", mp, " -mr ", mr, " -pc ", pc, " -pt ", pt, " -gold ", gold, " -ty SMUG", sep="")
#	system(command)
#	setwd(cur)
#}
