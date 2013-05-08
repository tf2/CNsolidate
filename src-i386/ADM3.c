/*
	An interpretation of the ADM2 - automated detection algorithm.	
	Author:	Tomas William Fitzgerald
	Date:	13/03/2011	
	
*/

#include <R.h>
#include <Rmath.h>
#include <math.h>

void _SADM(double *d, double *dat, double *prob, int *size);
void _LADM(int *start, int *stop, double *score, int *size);
void _EADM(double *d, double *dat, double *prob, int *size, double* t, double* rd);
//void _LADM(int *s, int *ss, double *r, double *w, int *n, int *start, int *stop, double *ratio, double *weight, double *score, int *size);

double sum3(double *array, int start, int stop) { double sum=0; for(int i=start;i<=stop;i++) sum+=array[i]; return sum; }
double msum3(double *array1, double *array2, int start, int stop) { double sum=0; for(int i=start;i<=stop;i++) { double a=array1[i]*array2[i]; sum+=a; } return sum; }
double mean3(double *array, int N) { if(N==1) { return array[0]; } else { double sum=0; for (int i=0; i<N; i++){ sum=sum+array[i]; } return sum/N; } }
double mmean3(double *array, int start, int stop) { int N=stop-start; if(N==0) { return array[0]; } else { double sum=0; for (int i=start; i<=stop; i++){ sum=sum+array[i]; } return sum/N; } }
double var3(double *array, int N) { int D = N-1; double xbar[N]; double m = mean3(array, N); double var = 0; for (int i=0;i<N;i++) { double dev = array[i] - m; var += pow(dev,2); } return(var/D); } 

void _SADM(double *d, double *dat, double *prob, int *size) {
	double v=var3(dat, *size); double m= mean3(dat, *size);
	for(int i=0;i<*size;i++) { dat[i] =(dat[i]-m)/v; prob[i] = 1/(pow(prob[i],2)); d[i] = (prob[i]*dat[i])/sqrt(prob[i]); }
}

void _LADM(int *start, int *stop, double *score, int *size) {
	int pin=0;
	for(int i=0;i<*size-1;i++) {
		start[pin]=i+1;
		if(score[i]==score[i+1]) { while(score[i]==score[i+1]) i++; stop[pin]=i+1; } else { stop[pin]=i+1; } pin++;	
	}	
}

void _EADM(double *d, double *dat, double *prob, int *size, double* t, double *rd) {	
	int i=0; int j=0;
	do {
		if(j>=*size | i >=*size) { break; }
		double q=sum3(prob,i, j); double s=msum3(prob, dat, i, j);
		double sc=s/sqrt(q); d[j]=sc;
		if(abs(sc)>*t) { 
			double scc=sc; int c=0;
				while(abs(sc)>*t & j<*size & abs(scc) <= abs(sc) & abs(dat[j]) >= abs(dat[j-1]/ *rd)) {
					j++; c++; scc=sc;
					if(j==*size | i ==*size) { break; }
					q=sum3(prob,i, j); s=msum3(prob, dat, i, j); sc=s/sqrt(q);  	
			 	}		
			 if(c==0) { i++; j++; } if(c==1) { i=j; }
			 if(c>1) { for(int k=i;k<j;k++) { d[k]=sc; } j++; i=j; }
		}
		else { i++; j++; }
	} while(j<*size & i <*size);
}

/* NB. I want to speed up the indexing by doing it in C... but its a bit of a pain - see the R call to this function to learn more.
void _LADM(int *s, int *ss, double *r, double *w, int *n, int *start, int *stop, double *ratio, double *weight, double *score, int *size) {
	int pin=0; int p1[*size]; int p2[*size];
	for(int i=0;i<*size-1;i++) {
		p1[pin]=i;
		if(score[i]==score[i+1]) { while(score[i]==score[i]) { i++; } p2[pin]=i; } else { p2[pin]=i; } pin++;	
	} int count = pin; pin=0;
	for(int i=0;i<count;i++) {
		if(p1[i]!=0 & p2[i]!=0) { s[pin]=start[p1[i]]; ss[pin]=stop[p2[i]]; r[pin]=mmean3(ratio, p1[i], p2[i]); w[pin]=mmean3(weight, p1[i], p2[i]);  n[pin]=(p2[i]-p1[i])+1;  pin++; }
	}
}
*/