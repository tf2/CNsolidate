######################################################
##             INITIALIZE METHOD                    ##
######################################################

## Initialize method for class SMAPProfile
setMethod("initialize", "SMAPProfile",
          function(.Object,
                   HMM,
                   observations,
                   P,
                   Q,
                   name=character(0)) {

              ## Store values in object
              .Object@HMM <- HMM
              .Object@observations <- observations
              .Object@P <- P
              .Object@Q <- Q
              .Object@name <- name

              ## Check lengths of required vectors
              len <- c(noObservations(observations),
                       length(Q))
              len.un <- unique(len)
              if (length(len.un) > 1)
                  stop("observations and state sequence of unequal lengths")
              if (len.un == 0)
                  stop("observations and state sequence must be of length > 0")

              ## Check state consistency
              state.un <- unique(Q)
              nas <- is.na(state.un)
              if (any(!state.un[!nas] %in% 1:noStates(HMM)))
                  stop("Unrecognized states in state sequence")

              ## Return new object
              .Object
          })

SMAPProfile <- function(HMM, observations, P, Q, name=character(0)) {

    new("SMAPProfile", HMM=HMM, observations=observations,
        P=P, Q=Q, name=name)
}

setMethod("show", "SMAPProfile", function(object) {

    cat("An object of class \"SMAPProfile\"\n")

    cat(paste("\nP:", P(object), "\n"))

    cat(paste("Use methods HMM(object), observations(object),",
              "P(object), and Q(object)",
              "to access object slots.\n"))
})

setMethod("show", "SMAPProfiles", function(object) {

    cat("An object of class \"SMAPProfiles\"\n")

    cat(paste("\nName:", name(object), "\n"))

    show(object@.Data)
})

######################################################
##             ACCESSOR METHODS                     ##
######################################################

setMethod("HMM", "SMAPProfile", function(object) object@HMM)
setMethod("observations", "SMAPProfile", function(object) object@observations)
setMethod("P", "SMAPProfile", function(object) object@P)
setMethod("Q", "SMAPProfile", function(object) object@Q)
setMethod("name", "SMAPProfile", function(object) object@name)
setMethod("name", "SMAPProfiles", function(object) object@name)

setMethod("observations", "SMAPProfiles",
          function(object) {

              obs <- observations(object[[1]])
              obs@name <- name(object)

              no.profiles <- length(object)

              if (no.profiles > 1) {
                  for (i in 2:no.profiles) {

                      p.obs <- observations(object[[i]])

                      obs@value <- c(obs@value, p.obs@value)
                      obs@chromosome <- c(obs@chromosome, p.obs@chromosome)
                      obs@startPosition <- c(obs@startPosition,
                                             p.obs@startPosition)
                      obs@endPosition <- c(obs@endPosition, p.obs@endPosition)
                      obs@reporterId <- c(obs@reporterId, p.obs@reporterId)
                  }
              }
              obs@noObservations <- length(obs@value)
              obs@chroms <- unique(obs@chromosome)
              obs@chrom.start <- match(obs@chroms, obs@chromosome)

              obs
          })

setMethod("Q", "SMAPProfiles",
          function(object) {

              Q <- Q(object[[1]])

              no.profiles <- length(object)

              if (no.profiles > 1)
                  for (i in 2:no.profiles)
                      Q <- c(Q, Q(object[[i]]))

              Q
          })

setMethod("[", "SMAPProfile",
          function(x, i, j, ..., drop) {
              x@observations <- x@observations[i]
              x@Q <- x@Q[i]
              x
          })

######################################################
##          REPLACEMENT METHODS                     ##
######################################################

setReplaceMethod("[[", "SMAPProfiles",
                 function(x, i, j, value) {
                     if (!is(value, "SMAPProfile"))
                         stop("value must be of class SMAPProfile")
                     x@.Data[[i]] <- value
                     x
                 })

