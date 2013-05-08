`rSegmentation` <-
function(data, reso, SorB = 0) {
	
	u = unique(data[,1])
	print(u)
 	
	A <- list()
	for( x in 1:length(u)) {
		
		if (SorB == 0) {
			if(u[x] < 23) {
				A[[x]] = rId1(data, u[x], reso)
			}
		}
		
		if (SorB == 1) {
			A[[x]] = rId2(data, u[x], reso)
			print(x)
		}
		if (SorB == 2) {
			A[[x]] = rId(data, u[x], length(data[data[,1]==u[x],1]))
			print(x)
		}
		
	}
	return(A)
}

