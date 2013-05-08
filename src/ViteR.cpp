#include<set>
#include<map>
#include<math.h>
#include<vector>
#include<sstream>
#include<fstream>
#include<iostream>

using namespace std;

#include<R.h>
#include<Rdefines.h>
#include<Rinternals.h>
#include<Rinterface.h>

extern "C" {

const double Pi = 3.141592;
map < string, vector<string> > mdata;

class hmm {
	public:
    	hmm(): d_nsts(0), d_nem(0), d_ntr(0), d_emsp(0), d_trsp(0) {}
    	hmm(double nsts, int nem, int ntr, double *emsp, double *trsp): d_nsts(nsts), d_nem(nem), d_ntr(ntr), d_emsp(emsp), d_trsp(trsp){}
	public:
   		int d_nsts; int d_nem; int d_ntr; double *d_emsp; double *d_trsp;
};

double npdf(double  mu, double sigma, double x) {
	double f = (x - mu)/sigma; double p = exp(-0.5*pow(f,2) ) / (sqrt(2*Pi) * sigma);
return p;
}

vector < double** > ets(hmm h) {
	double **ep; double **tp; 
	int N = (int)h.d_nsts; ep = new double*[N]; tp = new double*[N]; int c=0;
	for (int i=0;i<N;++i) { ep[i] = new double[2]; tp[i] = new double[3]; }
	for(int i=0;i<N;i++) {
		for(int j=0;j<h.d_nem;j++) { ep[i][j]=h.d_emsp[c]; c++; }
	} c=0;
	for(int i=0;i<N;i++) {
		for(int j=0;j<h.d_ntr;j++) { tp[i][j]=h.d_trsp[c]; c++; }
	}
	vector<  double** > r; r.push_back(ep); r.push_back(tp);
return r;
}

vector<double> viteR(hmm h, vector<double> d, vector<double> s) {        
	int i, j; int N = (int)h.d_nsts; 
	int **fwds = new int*[N]; int **bwds = new int*[N]; 
	vector < double** > rrr = ets(h); double **ep = rrr[0]; double **tp = rrr[1];
    for(i=0;i<N;++i) { fwds[i] = new int[d.size()]; bwds[i] = new int[d.size()]; }
    for(i=0;i<N;i++) { fwds[i][0] = log(npdf(ep[i][0], ep[i][1], d[0])); }  
        for (int i=1;i<d.size();i++) {
            for (int j=0;j<N;j++) {
                int pps; double mp; double sps[N];
                for (int ps = 0; ps < N; ps++) { sps[ps] = fwds[ps][i-1] + log(tp[ps][j]); }
                mp = sps[0];  pps = 0;
				for (int k=0;k<N;k++) {
					if (sps[k] > mp) { mp = sps[k]; pps = k; }
				}
                fwds[j][i] = mp + log(npdf(ep[j][0], ep[j][1], d[i])); bwds[j][i] = pps;
         	}
        }  
        int fs=0; double ml=fwds[0][d.size()-1];
        if(fwds[1][d.size()-1]>ml){ ml=fwds[1][d.size()-1]; fs=1; }
        if(fwds[2][d.size()-1]>ml){ ml=fwds[2][d.size()-1]; fs=2; }
        s[d.size()-1]=fs; int ns=fs;
        for(int i= d.size()-2;i>=0;i--){ s[i]=bwds[ns][i+1]; ns=s[i]; }
return s;        
}

vector<double> nord(double *data, int N) {
	double su=0; int count=0; vector < double > e;
	for(int i=0;i<N;i++) { su+=data[i]; count++; }
	double sd; double sus=0; double m=su/(count);    
    for (int i=0;i<count; i++) {  sus+=pow((data[i]-m),2); }
    double var = sus/(count-1); sd=sqrt(var);
    for (int i=0; i<count;i++) { data[i]=(data[i]- m)/(float)sd; e.push_back(data[i]); }
return e;
}

vector< double > maprecip(int a1, int a2, int b1, int b2) {
	double d; double dd; vector< double > recip;
	if(a2<b1 | a1>b2) { recip.push_back(0); recip.push_back(0); return recip; }
	if(a1==b1 & a2 == b2) { recip.push_back(1); recip.push_back(1); return recip; }
	if(a1<=b1 & a2>b2) { d=b2-b1; dd=a2-a1; d=(d/dd); recip.push_back(d); recip.push_back(1); return recip; }
	if(a1>=b1 & a2<b2) { recip.push_back(1); d=a2-a1; dd=b2-b1; d=(d/dd); recip.push_back(d); return recip; }
	if(a1==b2 & a2 > b2) { d=1; dd=a2-a1; d=(d/dd); recip.push_back(d); d=1; dd=b2-b1; d=(d/dd); recip.push_back(d); return recip; }
	if(a2==b1 & a1 < b1) { d=1; dd=a2-a1; d=(d/dd); recip.push_back(d); d=1; dd=b2-b1; d=(d/dd); recip.push_back(d); return recip; }
	if(a1<b1 & a2>b1 & a2<=b2) { d=min(a2-b1, a2-a1); dd=max(a2-b1, a2-a1); d=(d/dd); recip.push_back(d); d=min(a2-b1, b2-b1); dd=max(a2-b1, b2-b1); d=(d/dd); recip.push_back(d); return recip; }
	if(a2>b2 & a1<b2 & a1>=b1) { d=min(a2-a1, b2-a1);  dd=max(a2-a1, b2-a1); d=(d/dd); recip.push_back(d); d=min(b2-b1, b2-a1); dd=max(b2-b1, b2-a1); d=(d/dd); recip.push_back(d); return recip; }
	recip.push_back(0); recip.push_back(0); return recip;
}
	 
void mapfreq(int a, int b, int *c, int *d, double *fre1, double *fre2, double *fre3, double *r1, double *r2, double *rf1, double *rf2, double *rf3, double *ty, double *rty, int i, int N) {
	double r = 0; double rr=0; double f1=-1; double f2=-1; double f3=-1; double t1=ty[0];
	for(int j=0;j<N;j++) {
		vector < double > recip = maprecip(a,b,c[j],d[j]);
		if ( (recip[0]+ recip[1]) > (r+rr) ) {
			r=recip[0]; rr=recip[1]; f1=fre1[j]; f2=fre2[j]; f3=fre3[j]; t1=ty[j];
		}
	}
	r1[i] = r; r2[i] = rr;
	rf1[i] = f1; rf2[i] = f2; rf3[i] = f3; rty[i]=t1;
}
	 
void mapfrequency(int *a, int *b, int *c, int *d, double *fre1, double *fre2, double *fre3, double *r1, double *r2, double *rf1, double *rf2, double *rf3, double *ty, double *rty, int *size1, int *size2) {
	int N=(int)*size1; int NN=(int)*size2;
	for( int i=0;i<N;i++) {
		mapfreq((int)a[i], (int)b[i], c, d, fre1, fre2, fre3, r1, r2, rf1, rf2, rf3, ty, rty, i, NN);
	}
}
	
void ViteR(double *data, double *states, double *emissions, double *transitions, int *dN, int *sN, int *eN, int *tN) {
	vector < double > d = nord(data, (int)*dN);
	vector < double > s; for(int i=0;i<(int)*dN;i++) { s.push_back(states[i]); }
	hmm h((int)*sN, (int)*eN, (int)*tN, emissions, transitions);
	vector< double > r = viteR(h,d,s);
	for(int i=0;i<(int)*dN;i++) {
		states[i]=r[i];
	}
}

}
