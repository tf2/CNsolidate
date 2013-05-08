`selectHyper` <-
function(obs, p0, q0, K=20, M=10, iterMax=20) {
   param = .C("selectHyper", 
               obs, 
               as.integer(length(obs)), 
               as.double(p0), 
               as.double(q0), 
               as.integer(K),
               as.integer(M),
               as.integer(iterMax), 
               output = double(6),
               PACKAGE = "cnvSource")$output
    
   list(p=param[1], a=param[2], b=param[3], mu=param[4], 
          v=param[5], sigma2=param[6])
}

