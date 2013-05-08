`wscore` <- function(d) {
	wf = 0
	if (length(d) > 5000) {
		l = round(length(d)/10)
       	ll = round(length(d)/500)
       	if(!l%% 2) { l=l+1; }
       	if(!ll%% 2) { ll=ll+1; }
       	rem = runmed(d, ll)-runmed(d, l)
       	rem=rem[abs(rem)<quantile(abs(rem), probs=0.68)]
       	wf = (max(rem)-min(rem))*14
    }
return(wf)
}
