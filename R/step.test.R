`step.test` <- function(data=NULL) {
 	V <- function(f, len) {
		x <- seq (-1, 1, length = len); xx = x^(len*2)+(x^2)*f;	
		xs = ((xx- (1*min(x)))/((1*max(x))-(1*min(x))))
	return(xs/max(xs))
	}
	s = seq(0,5, by=0.01)
	ss = vector()
	for(x in 1:length(s)) {
		v = V(s[x], length(data))
		v = (v-mean(v))+mean(data)
		ss[x] = sum(abs(data-v))
	}
	step.like =1-(which.min(ss)/length(s))
return(step.like)
}