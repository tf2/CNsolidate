`syn.example` <- function() {
	chrs = rbind(cbind(1, 10000,0.001,0.2,0.5),
		cbind(2, 10000,0.001,0.2,0.5),
		cbind(3, 10000,0.001,0.2,0.5),
		cbind(3, 10000,0.001,0.2,0.5),
		cbind(4, 10000,0.001,0.2,0.5),
		cbind(5, 10000,0.001,0.2,0.5),
		cbind(6, 10000,0.001,0.2,0.5),
		cbind(7, 10000,0.001,0.2,0.5),
		cbind(8, 10000,0.001,0.2,0.5),
		cbind(8, 10000,0.001,0.2,0.5),
		cbind(9, 10000,0.001,0.2,0.5),
		cbind(10, 10000,0.001,0.2,0.5),
		cbind(11, 10000,0.001,0.2,0.5),
		cbind(12, 10000,0.001,0.2,0.5),
		cbind(13, 10000,0.001,0.2,0.5),
		cbind(14, 10000,0.001,0.2,0.5),
		cbind(15, 10000,0.001,0.2,0.5),
		cbind(16, 10000,0.001,0.2,0.5),
		cbind(17, 10000,0.001,0.2,0.5),
		cbind(18, 10000,0.001,0.2,0.5),
		cbind(19, 10000,0.001,0.2,0.5),
		cbind(20, 10000,0.001,0.2,0.5),
		cbind(21, 10000,0.001,0.2,0.5),
		cbind(22, 10000,0.001,0.2,0.5)
		)


m = syn.genome(chrs)
d = m$data
d[,1] = paste("chr", d[,1], sep="")
write.table(d, file="data.bed", sep="\t", row.names=F, col.names=F, quote=F)
write.table(m$rep, file="gold_report.txt", sep="\t", row.names=F, col.names=F, quote=F)

}