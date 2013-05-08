`smooth.CNA` <- function(x, smooth.region=2, outlier.SD.scale=4,
                       smooth.SD.scale=2, trim=0.025)
  {
    if (!inherits(x, 'CNA')) stop("First arg must be of class CNA")
    nsample <- ncol(x)-2
    chrom <- x$chrom
    uchrom <- unique(chrom)
    if(attr(x, "data.type")=="binary") stop("Not smoothing binary data ")
    for (isamp in 1:nsample) {
      genomdat <- x[,isamp+2]
      ina <- which(!is.na(genomdat) & !(abs(genomdat)==Inf))
      trimmed.SD <- sqrt(trimmed.variance(genomdat[ina], trim))
      outlier.SD <- outlier.SD.scale*trimmed.SD
      smooth.SD <- smooth.SD.scale*trimmed.SD
      k <- smooth.region
      for (i in uchrom) {
        ina <- which(!is.na(genomdat) & !(abs(genomdat)==Inf) & chrom==i)
        n <- length(genomdat[ina])
        smoothed.data <- .Fortran("smoothLR",
                                  as.integer(n),
                                  as.double(genomdat[ina]),
                                  sgdat=double(n),
                                  as.integer(k),
                                  as.double(outlier.SD),
                                  as.double(smooth.SD),
                                  PACKAGE = "CNsolidate")$sgdat
        x[,isamp+2][ina] <- smoothed.data
      }
    }
    x
  }

