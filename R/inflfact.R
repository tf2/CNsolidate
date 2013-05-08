`inflfact` <- function(trim)
  {
    a <- qnorm(1-trim)
    x <- seq(-a,a,length=10001)
    x1 <- (x[-10001] + x[-1])/2
    1/(sum(x1^2*dnorm(x1)/(1-2*trim))*(2*a/10000))
  }