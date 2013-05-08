#include <R.h>
#include <Rmath.h>

void _MAP(double *dat, int *ind1, int *ind2, double *rat, int *size);
double meanM(double *array, int start, int stop) { int N=stop-start; if(N==0) { return array[start]; } else { double sum=0; for(int i=start;i<=stop;i++) { sum=sum+array[i]; } return sum/N; } }

void _MAP(double *dat, int *i1, int *i2, double *rat, int *size) {
	int N=(int)*size; double mr = (double)*rat; double m= meanM(dat, i1[0], i2[0]);	
	
	if(m>0) {
		//if(i1[0]<N & i1[0]>0) { while(dat[i1[0]]<mr) { i1[0]++; if(i1[0]>N | i1[0]<1) { break; } } }
		if(i1[0]<N & i1[0]>0) { while(dat[i1[0]]<mr | dat[i1[0]]<m/2) { i1[0]++; if(i1[0]>N | i1[0]<1) { break; } } }
		if(i1[0]<N & i1[0]>1) { while(dat[i1[0]-1]>mr & dat[i1[0]-1]>m/2) { i1[0]--; if(i1[0]>=N | i1[0]<1) { break; } } }
		//if(i1[0]<N & i1[0]>1) { while(dat[i1[0]-1]>mr & dat[i1[0]-1]>atan(m)/1.62) { i1[0]--; if(i1[0]>=N | i1[0]<1) { break; } } }
		if(i2[0]<N & i2[0]>1) { while(dat[i2[0]]<mr) { i2[0]--; if(i2[0]>N | i2[0]<1) { break; } } }
		if(i2[0]<N & i2[0]>0) { while(dat[i2[0]+1]>mr & dat[i2[0]+1]>m/2) { i2[0]++; if(i2[0]>=N | i2[0]<1) { break; } } }
		//if(i2[0]<N & i2[0]>0) { while(dat[i2[0]+1]>mr & dat[i2[0]+1]>atan(m)/1.62) { i2[0]++; if(i2[0]>=N | i2[0]<1) { break; } } }
	}
	if(m<0) {
		//if(i1[0]<N & i1[0]>0) { while(dat[i1[0]]>-mr) { i1[0]++; if(i1[0]>=N | i1[0]<1) { break; } } }
		if(i1[0]<N & i1[0]>0) { while(dat[i1[0]]>-mr | dat[i1[0]]>m/2) { i1[0]++; if(i1[0]>=N | i1[0]<1) { break; } } }
		if(i1[0]<N & i1[0]>1) { while(dat[i1[0]-1]<-mr & dat[i1[0]-1]<m/2) { i1[0]--; if(i1[0]>=N | i1[0]<1) { break; } } }
		//if(i1[0]<N & i1[0]>1) { while(dat[i1[0]-1]<-mr & dat[i1[0]-1]<atan(m)/1.62) { i1[0]--; if(i1[0]>=N | i1[0]<1) { break; } } }
		if(i2[0]<N & i2[0]>1) { while(dat[i2[0]]>-mr) { i2[0]--; if(i2[0]>N | i2[0]<1) { break; } } }
		if(i2[0]<N & i2[0]>0) { while(dat[i2[0]+1]<-mr & dat[i2[0]+1]<m/2) { i2[0]++; if(i2[0]>=N | i2[0]<1) { break; } } }
		//if(i2[0]<N & i2[0]>0) { while(dat[i2[0]+1]<-mr & dat[i2[0]+1]<atan(m)/1.62) { i2[0]++; if(i2[0]>=N | i2[0]<1) { break; } } }
	} rat[0]= meanM(dat, i1[0], i2[0]);
}
