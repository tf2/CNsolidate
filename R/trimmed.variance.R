`trimmed.variance` <- function(genomdat, trim=0.025)
  {
    n <- length(genomdat)
    n.keep <- round((1-2*trim)*(n-1))
    inflfact(trim)*sum((sort(abs(diff(genomdat)))[1:n.keep])^2 / (2*n.keep))
  }