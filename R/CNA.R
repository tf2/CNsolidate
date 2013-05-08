`CNA` <- function(genomdat, chrom, maploc, data.type=c("logratio","binary"),
                   sampleid=NULL)
  {
    if (!is.numeric(genomdat)) stop("genomdat must be numeric")
    if (is.factor(chrom)) chrom <- as.character(chrom)
    if (!is.numeric(maploc)) stop("maploc must be numeric")
    data.type <- match.arg(data.type)
    ina <- (!is.na(chrom) & is.finite(maploc))
    if (sum(!ina)>0)
      warning("markers with missing chrom and/or maploc removed\n")
    sortindex <- which(ina)[order(chrom[ina], maploc[ina])]
    if (is.vector(genomdat)) genomdat <- as.matrix(genomdat)
    if (!missing(sampleid)) {
      if (length(sampleid) != ncol(genomdat)) {
        warning("length(sampleid) and ncol(genomdat) differ, names ignored\n")
        sampleid <- paste("Sample", 1:ncol(genomdat))
      } 
    } else {
        sampleid <- paste("Sample", 1:ncol(genomdat))
    }
    colnames(genomdat) <- sampleid
    zzz <- data.frame(chrom=I(chrom), maploc=maploc, genomdat)
    zzz <- zzz[sortindex,]

    if(!all(sapply(unique(chrom),function(ichrom,chrom,maploc){
      length(maploc[chrom==ichrom])-length(unique(maploc[chrom==ichrom]))
    },chrom,maploc)==0)) warning("array has repeated maploc positions\n")

    attr(zzz, "data.type") <- data.type
    class(zzz) <- c("CNA","data.frame")
    zzz
  }





  
  
  

  