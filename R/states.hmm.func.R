states.hmm.func <-
    function(sample, chrom, dat, datainfo = clones.info, vr = .01,
             maxiter = 100, aic = FALSE, bic = TRUE, delta = 1,
             nlists = 1, eps = .01, print.info = FALSE,
             diag.prob = .99)
{

    obs <- dat[datainfo$Chrom==chrom, sample]
    kb <- datainfo$kb[datainfo$Chrom==chrom]
    ##with current sproc files, data is already ordered by kbs
    obs.ord <- obs[order(kb)]
    kb.ord <- kb[order(kb)]

    ind.nonna <- which(!is.na(obs.ord))

    y <- obs.ord[ind.nonna]
    kb <- kb.ord[ind.nonna]


#####################################

    numobs <- length(y)
    zz <- vector(mode = "list", 5)
    zz[[1]] <-
        list(log.lik =
             sum(dnorm(y, mean = mean(y), sd = sd(y), log = TRUE)))
    for(k in 2:5)
    {
        
        mu <- kmeans(y, k)$centers
        gamma <- matrix((1 - diag.prob) / (k - 1), k, k)
        diag(gamma) <- diag.prob
        zz[[k]] <-
        {

            res <-
                .C("calc_observed_likelihood",
                   as.integer(numobs),
                   as.double(y),
                   as.integer(k),
                   mu = as.double(mu),
                   sigma = as.double(sqrt(vr)),
                   gamma = as.double(gamma),
                   pi = as.double(rep(-log(k), k)),
                   num.iter = as.integer(maxiter),
                   as.double(eps),
                   log.lik = double(1),
                   filtered.cond.probs = double(k * numobs),
                   hidden.states = integer(numobs),
                   as.logical(print.info),
                   PACKAGE = "CNsolidate")
            res$hidden.states <- res$hidden.states + 1
            res$filtered.cond.probs <-
                matrix(res$filtered.cond.probs, nr = k)
            res$gamma <- matrix(res$gamma, nr = k)
            res
            
        }
        
    }

###############################################3
###############################################3
    ##identify the model with the smallest model selection criteria

    ##now, scroll over all options:

    ##number of states (means) + number of states*(number of states-1) (transitions) #+ 1 (variance)

#    kk <- c(2, 5, 10, 17, 26)
    kk <- (1:5) ^ 2 + 1
    for (nl in 1:nlists)
    {
        if ((aic) && (nl==1))
        {
            ##-loglik+2*k/2
            factor <- 2
        }
        else if (bic)
        {
            ##-loglik+log(n)*k*delta/2
            if (aic)
            {
                factor <- log(numobs)*delta[nl-1]
            }
            else
            {
                factor <- log(numobs)*delta[nl]
            }
        }
        lik <- sapply(zz, function(z) -z$log.lik) + kk * factor / 2
        nstates <- likmin <- which.min(lik)
        z <- zz[[likmin]]

######################################
        ##out rpred and state

        if (nstates > 1) #if non-generic
        {
            ##print(nstates)
            maxstate <- apply(z$filter, 2, which.max)
###            maxstate <- z$hidden.states
            rpred <- as.vector(z$mu %*% z$filter)
            prob <- apply(z$filter, 2, max)
            ##use median for prediction and mad for state dispersions
            maxstate.unique <- unique(maxstate)
            pred <- rep(0, length(y))
            disp <- rep(0, length(y))
            for (m in 1:length(maxstate.unique))
            {
                
                pred[maxstate==maxstate.unique[m]] <-
                    median(y[maxstate==maxstate.unique[m]])
                disp[maxstate==maxstate.unique[m]] <-
                    mad(y[maxstate==maxstate.unique[m]])
                
            }

            ##if (length(z$pshape) == 1)
            ##{
            ##        disp <- rep(z$pshape, length(maxstate))
            ##}
            ##else
            ##{
            ##        disp <- z$pshape[maxstate]
            ##}
            
        }
        else #if generic
        {
            maxstate <- rep(1, length(y))
            ##rpred <- rep(mean(y), length(y))
            rpred <- rep(median(y), length(y))
            prob <- rep(1, length(y))
            ##pred <- rep(mean(y), length(y))
            pred <- rep(median(y), length(y))
            ##disp <- rep(var(y), length(y))
            disp <- rep(mad(y), length(y))
            
        }
        
        out <-
            cbind(matrix(maxstate, ncol=1),
                  matrix(rpred, ncol=1),
                  matrix(prob, ncol=1),
                  matrix(pred, ncol=1),
                  matrix(disp, ncol=1))
        
        out.all <- matrix(NA, nrow=length(kb.ord), ncol=6)
        out.all[ind.nonna,1:5] <- out
        
        out.all[,6] <- obs.ord
        out.all <- as.data.frame(out.all)
        dimnames(out.all)[[2]] <- c("state", "rpred", "prob", "pred", "disp", "obs")
        
        
        if (nl==1)
        {
            out.all.list <- list(out.all)
            nstates.list <- list(nstates)
        }
        else
        {
            out.all.list[[nl]] <- out.all
            nstates.list[[nl]] <- nstates
        }
        
        ##cloneinfo <- as.data.frame(cbind(rep(chrom, length(kb.ord)), kb.ord))
        ##dimnames(cloneinfo)[[2]] <- c("Chrom", "kb")
    }
    list(out.list = out.all.list, nstates.list = nstates.list)
    
}
