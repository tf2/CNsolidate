`getBcmixSmoothClass` <-
function(obs, hyper, pos=1:length(obs), 
    classify="F", K=20, M=10, thres=sqrt(as.double(hyper[6]))*2) {

    fileNameRand = round(runif(1)*1000000)
    fit = .C("cppBcmixSmooth",
       obs,
       as.integer(length(obs)),
       as.double(hyper[1]),
       as.double(hyper[2]),
       as.double(hyper[3]),
       as.double(hyper[4]),
       as.double(hyper[5]),
       as.double(hyper[6]),
       as.integer(ifelse(classify=="F",0,1)),
       as.integer(K),
       as.integer(M),
       estSig = double(length(obs)),
       estNon0prob = double(length(obs)),
       PACKAGE = "CNsolidate");
      
       sig = fit$estSig;
       non0prob = fit$estNon0prob;


       stateseq=NULL;
       if (classify == "T") { 
            # Classify each location as "G", "L", or "0" and store in stateseq.

       } 

       zzz = list(obs=obs, pos=pos, hyper=hyper, sig=sig, non0prob=non0prob, stateseq=stateseq, K=K, M=M)
       class(zzz) = c("CNsolidate")
       zzz
}

