`smugFun` <-
function(l) {
	sm = list(p=0.95, a=0.99, b=0.005, mu=0, v=1, sigma2=var(l))
	obs.bcmix<-getBcmixSmoothClass(l, hyper=sm, classify="T", thres=0.1)
	return(obs.bcmix$non0prob)
}

