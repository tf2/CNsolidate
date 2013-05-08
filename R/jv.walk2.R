`jv.walk2` <- function (file, idir, odir, reso = 1000, pc = 0.05, p1 = 0.95, 
    p2 = 0.01, mp = 5, mr = 0.2, pt = 3, gold = "../data/42Mcalls_all_feb09.txt") 
{
    cur = getwd()
    Jpath = system.file("java", package = "CNsolidate")
    setwd(Jpath)
    d = read.table(file)
    n = unlist(strsplit(file, "/"))
   	f =  paste(odir, substr(n[length(n)], 0, nchar(n[length(n)]) - 4), "_variancewalking.temp", sep = "")
   	name = paste(odir, substr(n[length(n)], 0, nchar(n[length(n)]) - 4), "_variancewalking_walk.temp", sep = "")
    rep = paste(odir, substr(n[length(n)], 0, nchar(n[length(n)]) - 4), "_variancewalking_FinalReport_walk.txt", sep = "")
    pfrep = paste(odir, "P",substr(n[length(n)], 0, nchar(n[length(n)]) - 4), "_FinalReport_walk.txt", sep = "")
    nfrep = paste(odir, "N",substr(n[length(n)], 0, nchar(n[length(n)]) - 4), "_FinalReport_walk.txt", sep = "")
    u = unique(d[,1])
    v = vector()
    x=1
    dd=NULL
    for(x in 1:length(u)) {
    	if(x < 23) {	
    		write.table(d[d[,1]==u[x],], file=f, row.names=F, col.names=F, quote=F, sep="\t")
    		command = paste("java -Xmx1600m ControlT -f ", f, " -mr ", mr, " -i ", idir, " -o ", odir, sep = "")
    		system(command)
    		command = paste("java -Xmx1600m Clean -f ", name, " -mp ", mp, " -mr ", mr, " -pc ", pc, " -pt ", pt, " -gold ", gold, " -ty walk", sep = "")
    		system(command) 
    		if(length(count.fields(rep))>0) {
    			dd=rbind(dd, read.table(rep))
    		}
    	}
    }
    setwd(cur)
    command = paste("rm", f, rep, name, sep=" ")
    system(command)
    write.table(dd[dd[,4]>0,], file=pfrep, row.names=F, col.names=F, sep="\t", quote=F)
    write.table(dd[dd[,4]<0,], file=nfrep, row.names=F, col.names=F, sep="\t", quote=F)
}