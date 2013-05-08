`ADM_feature` <- function(ddd, np) {
      rrr <- .C("_LADM"
      			, start = as.integer(vector(length = length(ddd[,1])))
      			, stop = as.integer(vector(length = length(ddd[,1])))
      			, score = as.double(ddd[, 8])
      			, size = as.integer(length(ddd[,1]))
      			, PACKAGE = "CNsolidate")
      			
    c = cbind(rrr$start, rrr$stop); c = c[c[, 1] != 0 & c[, 2] != 0, ]; c = c[abs(c[, 2] - c[, 1]) + 1 >= np, ]
    v = vector(); vv = vector(); vvv = vector(); vvvv = vector();
    if (!is.null(dim(c)) & length(c)!=0) {
        for (x in 1:length(c[,1])) {
        	if(c[x,2]>length(ddd[,1])) { c[x,2] = length(ddd[,1]); }
            v[x] = mean(ddd[c[x,1]:c[x,2],4])
            vv[x] = mean(ddd[c[x,1]:c[x,2],7])
            vvv[x] = mean(ddd[c[x,1]:c[x,2],8])
            vvvv[x] = length(ddd[c[x,1]:c[x,2],8])
        }
    return(cbind(ddd[c[,1],1], ddd[c[,1],2], ddd[c[,2],3], v, vv, vvv, vvvv, c))
    }
    if(is.null(dim(c)) & length(c)==2) {
    	if(c[2]>length(ddd[,1])) { c[2] = length(ddd[,1]); }
    	v[1] = mean(ddd[c[1]:c[2],4])
        vv[1] = mean(ddd[c[1]:c[2],7])
        vvv[1] = mean(ddd[c[1]:c[2],8])
        vvvv[1] = length(ddd[c[1]:c[2],8])       
	return(cbind(ddd[c[1],1], ddd[c[1],2], ddd[c[2],3], v, vv, vvv, vvvv, c[1], c[2]))
	}
}