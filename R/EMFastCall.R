`EMFastCall`<-function(mdata,thr0,Factor)

{



        
        StartPar <- StartCond1(mdata, thr0)
    muvec <- StartPar$muvec
    sdvec <- StartPar$sdvec
    prior <- c(0.05, 0.1, 0.7, 0.1, 0.05)
    lvec <- c(-16, -0.9, -thr0, thr0, 0.7)
    uvec <- c(-0.9, -thr0, thr0, 0.7, 16)
    g = PosteriorP(mdata, muvec, sdvec, prior)
    g[!complete.cases(g),] = 0
    LikeliNew <- sum(g *  prior)
    threshold <- 1e-05
    
    for (i in 1:1000) {
        muvecold <- muvec
        LikeliOld <- LikeliNew
        taux <- EStep(mdata, muvec, sdvec, prior, lvec, uvec)
        MResult <- MStep(mdata, taux, muvec, sdvec)
        muvec <- MResult$mu
        sdvec <- MResult$sdev
        prior <- MResult$prior
        truncResult <- truncUp(muvec, sdvec, thr0, Factor)
        lvec <- truncResult$l
        uvec <- truncResult$u
        g = PosteriorP(mdata, muvec, sdvec, prior)
    	g[!complete.cases(g),] = 0
        LikeliNew <- sum(g * prior)
        if (abs(LikeliNew - LikeliOld) < threshold) {
            break
        }
    }
    
Results<-list()
Results$muvec<-muvec
Results$sdvec<-sdvec
Results$prior<-prior
Results$iter<-i
Results$bound<-cbind(lvec,uvec)
Results
}




