hmm.run.func <- function(data, vr = .01, maxiter = 100, aic = TRUE, bic = TRUE, delta = NA, eps = 0.01) {

	datainfo <- cbind(data[,1], data[,5] - data[,4])
	datainfo <- as.data.frame(datainfo)
	colnames(datainfo) <- c("Chrom", "kb")
	dat <- as.data.frame(data[,6])
	
    chrom.uniq <- unique(datainfo$Chrom)
    states <- matrix(NA, nrow=nrow(dat), ncol=(2+6*ncol(dat)))
    states[,1:2] <- cbind(datainfo$Chrom, datainfo$kb)
    nstates <- matrix(NA, nrow=length(chrom.uniq), ncol=ncol(dat))

    states.list <- list(states)
    nstates.list <- list(nstates)

    nlists <- 0
    if (aic){
        nlists <- 1
    }
    if (bic){
        if (is.na(delta)){
            delta <- c(1)
        }
        else{
            delta <- c(1, delta)
        }
        for (j in 1:length(delta)){
            nlists <- nlists+1
        }
    }
    
    if (nlists > 1){
        for (j in 2:nlists){
            states.list[[j]] <- states.list[[1]]
            nstates.list[[j]] <- nstates.list[[1]]
        }
    }
    
    for (i in 1:ncol(dat)) {
        
        colstart <- 2+(i-1)*6+1
        colend <- 2+i*6
        
        for (j in 1:length(chrom.uniq)){
                   
            res <- try(states.hmm.func(sample=i, chrom=j, dat=dat,
                                    datainfo=datainfo,vr=vr,
                                    maxiter=maxiter, aic=aic, bic=bic,
                                    delta=delta, nlists=nlists, eps = eps))
            if (!(inherits(res, "try-error")))
                for (m in 1:nlists) {

                    states.list[[m]][((1:nrow(states))[states[,1]==j]),colstart:colend] <-
                        as.matrix(res$out.list[[m]])
                    nstates.list[[m]][j,i] <- res$nstates.list[[m]]

                }
            else
                stop("hmm fit failed!\n")

        }
        
    }
    return(list(states.hmm = states.list, nstates.hmm = nstates.list))
}

