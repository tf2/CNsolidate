`rId` <- 
function (d, chr, reso) 
{
    t = 0
    chr1 = d[d[, 1] == chr, ]
    len = length(chr1[, 1])
    fac = ceiling(len/reso)
    st = 1
    sto = reso
    for (a in 1:fac) {
        dat = chr1[st:sto, 4]
        loc = chr1[st:sto, 1:3]
        dat = as.numeric(dat)
        if (length(dat) < 100) {
            zer = vector(length = 100)
            zzz = matrix(ncol=3, nrow=100)
            for (z in 1:length(zer)) {
                zer[z] = 0
                zzz[z,1] = loc[1,1]
                zzz[z,2] = loc[length(loc[,1]),2]+z
                zzz[z,3] = loc[length(loc[,1]),3]+z
            }
            temp = c(dat, zer)
            loc = rbind(loc,zzz)
            te = smapFun(temp, loc)
            pin = length(te) - 100
            t = c(t, te[1:pin])
        }
        else if (length(chr1[,1]) < reso) {
        	dat = chr1[1:len, 4]
        	t = c(t, smugFun(dat))
        }
        else {
            t = c(t, smapFun(dat, loc))
        }
        st = st + reso
        sto = sto + reso
        if (sto > len) {
            sto = len
        }
    }
    t = t[-(1)]
    return(cbind(chr1[, 1:4], t))
}

