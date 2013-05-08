`StartCond1`<-function(mdata,thr0)

{
muvec<- c(-1.5,-0.8,0,0.58,1)
sdvec<-c(0.01,0.01,0.01,0.01,0.01)
lvec<-c(-16,-0.9,-thr0,thr0,0.7)
uvec<-c(-0.9,-thr0,thr0,0.7,16)

for (i in 1:5)
{
u<-uvec[i]
l<-lvec[i]
ind<-which(mdata<=u & mdata>=l)
if (length(ind)==0)
{
muvec[i]<-muvec[i]
sdvec[i]<-sdvec[i]
}
if (length(ind)==1)
{
muvec[i]<-mdata[ind]
sdvec[i]<-sdvec[i]
}
if (length(ind)>1)
{
muvec[i]<-mean(mdata[ind])
sdvec[i]<-sd(mdata[ind])
}
}
Result<-list()
Result$muvec<-muvec
Result$sdvec<-sdvec
Result
}



