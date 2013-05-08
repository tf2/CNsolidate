`changepoints.prune` <- function(genomdat, lseg, change.cutoff=0.05) {
  n <- length(genomdat)
  nseg <- length(lseg)
  ncpt <- nseg-1
  zzz <- .Fortran("prune",
                  as.integer(n),
                  as.double(genomdat),
                  as.integer(nseg),
                  as.integer(lseg),
                  as.double(change.cutoff),
                  double(nseg),
                  as.integer(ncpt),
                  loc=integer(ncpt),
                  integer(2*ncpt),
                  pncpt=integer(1), PACKAGE="CNsolidate")
  pruned.ncpt <- zzz$pncpt
  pruned.cpts <- cumsum(lseg)[zzz$loc[1:pruned.ncpt]]
  pruned.lseg <- diff(c(0,pruned.cpts,n))
  pruned.lseg
}


