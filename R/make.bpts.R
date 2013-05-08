`make.bpts`<-function(x,pos,res,index){
	loc = which(x[,index]==pos)
	p = loc-res
	pp = loc+res
	l = length(x[,1])-res

	if(loc<res) {
		c =vector()
		for(i in 1:res) { c[i] = res-i }
		r1 = rep(x[loc,index],res)-c
		breaks = c(r1, x[1:pp,index])	
	} 
	
	if(loc>length(x[,1])-res) {
		c =vector()
		for(i in 1:res) { c[i] = res+i }
		r2 = rep(x[loc,index],res)+c
		breaks = c(x[l:length(x[,1]),index],r2)
	}
	
	if(!loc<res & !loc>l) {
		breaks<-x[p:pp,index]
	}
	#breaks<-breaks-1

return(sort(breaks))
}
