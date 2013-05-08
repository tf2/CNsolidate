`mergeAll` <-
function(lis, t=0.1, p1 =0.90, p2=0.1) {

	A <- list()
	for( x in 1:length(lis)) {
		#print(x)
		temp = merge1(lis[[x]], p=p1, pp=p2)
		temp[abs(temp[,6])<t,6] = 0
		A[[x]] = temp
	}
	
	return(A)
}

