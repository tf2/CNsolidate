`runAll` <- function(d, type=1, reso=1000, pc=0.1, p1 = 0.95, p2 = 0.01, mp=5, mr=0.3, pt = 3, gold = "42Mcalls_all_feb09.txt", out=d) {

	files = dir(d)
	
	for( x in 1:length(files)) {
		
		file = files[x]

		dt.g.wavelet(file, d, out)
		
		#safe.dec(file, d, out, reso, pc, p1, p2, mp, mr, pt, gold)
		#jv.walk(file, d, out, reso, pc, p1, p2, mp, mr, pt, gold)
		#SMUG(file, d, out, reso, pc, p1, p2, mp, mr, pt, gold)
		#GADA(file, d, out, pc, pt, mp, mr)
		#dna.copy(file, d, out, mp, mr)	
		#SMAP(file, d, out, reso, pc, p1, p2, mp, mr, pt, gold)
		#Fsegment(file, d, out, 1000, 2, mp, mr)
		#rFASTCALL(nfile, d, out)
		#runBCP(file, d, out, reso, pc, p1, p2, mp, mr, pt, gold)	
		
		pn.signs(file, d, out)
		
		pfile = paste("P", file, sep="")
		safe.dec(pfile, d, out, reso, pc, p1, p2, mp, mr, pt, gold)
		jv.walk(pfile, d, out, reso, pc, p1, p2, mp, mr, pt, gold)
		SMUG(pfile, d, out, reso, pc, p1, p2, mp, mr, pt, gold)
		GADA(pfile, d, out, pc, pt, mp, mr)
		dna.copy(pfile, d, out, mp, mr)	
		SMAP(pfile, d, out, reso, pc, p1, p2, mp, mr, pt, gold)
		Fsegment(pfile, d, out, 1000, 2, mp, mr)
		rFASTCALL(pfile, d, out)
		#runBCP(pfile, d, out, reso, pc, p1, p2, mp, mr, pt, gold)	
		
		#consen(pfile, d, out, 1, gold)

		nfile = paste("N", file, sep="")
		safe.dec(nfile, d, out, reso, pc, p1, p2, mp, mr, pt, gold)
		jv.walk(nfile, d, out, reso, pc, p1, p2, mp, mr, pt, gold)
		SMUG(nfile, d, out, reso, pc, p1, p2, mp, mr, pt, gold)
		GADA(nfile, d, out, pc, pt, mp, mr)
		dna.copy(nfile, d, out, mp, mr)		
		SMAP(nfile, d, out, reso, pc, p1, p2, mp, mr, pt, gold)
		Fsegment(nfile, d, out, 1000, 2,mp, mr)
		rFASTCALL(nfile, d, out)
		#rBCP(nfile, d, out, reso, pc, p1, p2, mp, mr, pt, gold)
		
		#consen(nfile, d, out, 1, gold)		
		#meeting.c.o(file, d, out)
		
		makeFObs(file, d, out, methods = c("safe", "walk", "SMUG", "GADA", "SMAP", "dna.copy", "fsegment", "fast"))
	}

}

