`gfct` <-function(x,moy,sdev,l,u)
{
	(dnorm(x,mean=moy,sd=sdev)*(x<=u)*(x>=l))/(pnorm(u,mean=moy,sd=sdev)-pnorm(l,mean=moy,sd=sdev))
}


