.onLoad <- function(lib, pkg) require("methods", quietly=TRUE)

## Class SMAPObservations
setClass("SMAPObservations",
         representation(## Vector of values (ratios or log2-ratios)
                        value="numeric",
                        ## Vector of chromosomes
                        chromosome="character",
                        ## Vector of unique chromosomes
                        chroms="character",
                        ## Vector of chromosome start positions
                        chrom.start="numeric",
                        ## Vector of start positions
                        startPosition="numeric",
                        ## Vector of end positions
                        endPosition="numeric",
                        ## Optional slots:
                        ## Assay name
                        name="character",
                        ## Vector of reporter ids
                        reporterId="character",
                        ## Derived slots
                        ## Vector of distance between clones
                        distance = "numeric",
                        ## Vector overlapping spots
                        overlapIds="numeric",
                        ## Vector of overlaps
                        overlaps="numeric",
                        ## Vector of start indices for each clone
                        ## in the previous two vectors
                        startOverlaps="numeric",
                        ## Vector of number of overlaps for each clone
                        noOverlaps="numeric",
                        noObservations="numeric"
                        ))

## Class SMAPHMM
setClass("SMAPHMM",
         representation(## Matrix of transition probabilities
                        A="matrix",
                        ## Vector of initial probabilities
                        Pi="numeric",
                        ## List of dparam objects for each distribution
                        Phi="list",
                        ## Length of states vector
                        noStates="numeric",
                        ## Matrix of transition probabilities
                        Z="matrix",
                        ## Vector of initial probabilities
                        Y="numeric",
                        eta="ANY",
                        grad="ANY"
                        ))

## Class SMAPProfile
setClass("SMAPProfile",
         representation(## SMAPHMM
                        HMM="SMAPHMM",
                        ## SMAPObservations
                        observations="SMAPObservations",
                        ## Joint posterior log probability
                        P="numeric",
                        ## State sequence
                        Q="numeric",
                        ## Name
                        name="character"
                        ))

## Class SMAPProfiles
setClass("SMAPProfiles",
         representation(name="character"),
         contains="list")

## Class grad
setClass("grad",
         representation(## Matrix of transition probabilities
                        A="matrix",
                        ## Vector of initial probabilities
                        Pi="numeric",
                        ## List of dparam objects for each distribution
                        Phi="list"
                        ))

## Class eta
setClass("eta",
         representation(value="numeric",
                        ## Matrix of transition probabilities
                        A="matrix",
                        ## Vector of initial probabilities
                        Pi="numeric",
                        ## List of dparam objects for each distribution
                        Phi="list"
                        ))

## Class gaussparam
setClass("gaussparam",
         representation(mean="numeric",
                        sd="numeric"))



validSegmentation = function(object) {
	
## check that slots 'breakpoints', 'logLik', 'hasConfint' all have the same length
	if(!((length(object@breakpoints)==length(object@logLik)) &&
		 (length(object@breakpoints)==length(object@hasConfint))))
    return(FALSE)
	
## check that nrow(y)==length(x)==length(flag)
	n = nrow(object@y)
	if(!(length(object@x)%in%c(0,n)))
    return(FALSE)
	if(!(length(object@flag)%in%c(0,n)))
    return(FALSE)
	
	if(!is.na(object@nrSegments))
    if(object@nrSegments<1 && object@nrSegments>length(object@breakpoints))
	return(FALSE)
	
	isGood = TRUE
## check the elements of the breakpoints slot (which is a list)
	for(i in seq(along=object@breakpoints)) {
		b = object@breakpoints[[i]]
		if(!((nrow(b)==i-1) && (ncol(b) %in% c(1,3)) &&
			 ("estimate" %in% colnames(b)) && all(b[, "estimate"] <= n, na.rm=TRUE)))
		isGood = FALSE
	}
	return(isGood)
}


setClass("segmentation",
representation(
y = "matrix",
x = "numeric",
flag = "integer",
breakpoints = "list",
logLik = "numeric",
hasConfint = "logical",
nrSegments = "integer"
),
prototype = list(
y = matrix(0, nrow=0, ncol=0),
x = numeric(0),
flag = integer(0),
breakpoints = list(),
logLik = numeric(0),
hasConfint = logical(0),
nrSegments = as.integer(NA)
),
validity = validSegmentation) ## see above


