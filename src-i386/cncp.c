#include <R.h>
#include <Rmath.h>
#include <math.h>

struct Boo {int a; int b; int c; int d; int bo;};

void wcncp(double *dat, double *prob, int *need, int *allow, int *index, int *index1, int *index2, int *size) {
	
	double pro = *prob; int p = *index; int n = *need; int al = *allow; int en = *size;
	struct Boo z; z.a=0; z.b=0; z.c=0; z.d=0; z.bo=0;
	
	while(z.bo==0) {
		
		if(dat[0]>pro) { z.a++; z.c=1; } 
		if(dat[0]==0 & z.c==1) { z.b++; } else { z.b=0; }
		if(z.a==n & z.d==0) { int t = p-1; *index1 = t; z.d=1; }
		
		if(z.b==al) { *index2 = p-1; z.bo=1; }
		if(p==en-1) { *index2 = p-1; z.bo=1; }
		dat++; p++;
	}
}

int cdoubles (const void *a, const void *b) {
    const double *da = (const double *) a;
    const double *db = (const double *) b;     
	return (*da > *db) - (*da < *db);
}

double quantile(double values[], double quantile) {
	double copy[sizeof(values)];
    for(int i=0;i<sizeof(copy);i++) {
        copy[i] =  fabs(values[i]);
    }
    qsort (copy, sizeof(copy), sizeof(double), cdoubles);
    int index = (int) (sizeof(copy) * quantile);
    return copy[index];
}

double *diff(double data[], int size) {
	double *d= malloc(sizeof(double)*size);
	for(int i=0;i<size-1;i++) {
		double ab = data[i]-data[i+1];
		if(ab>0) { d[i]=ab; } else { d[i]=ab/-1; }
		//  printf ("%f %f %f\n", data[i], data[i+1], d[i]);
	}
	d[size-1]=0;
	return(d);
}

double *mSum(double data[], int win, double fac) {
	
	int st = 0;
	int w = win;
	int nlen = sizeof(data)/win;
	double *cha = malloc(sizeof(double)*(win*nlen)+1);//[(win*nlen)+1];
	double val = quantile(data, 0.68)*fac;
	
	for(int i=0;i<nlen+1;i++) {
		double seg[w];
		int pin=0;
		for(int j=st;j<w;j++) {
			seg[pin]=data[j];
			pin++;
		}
		pin=0;
		double *tseg = diff(seg,win);
		for(int k=st;k<w;k++) {
			cha[k] = tseg[pin];
		}
		st = w+1;
		w=win*i;
	}
	
	for(int i=0;i<sizeof(cha);i++) {
		if(cha[i]>=val) { cha[i]=1; } else { cha[i]=0; }
	}

	return(cha);
}

double *cSum (double d[], int wind[], double fac, int size) {
	double *vec = malloc(sizeof(double)*size);//[size];
	double mat[size][sizeof(wind)/sizeof(int)];
	for(int x= 0;x<sizeof(wind);x++) {

		int c = (ceil(sizeof(d)/wind[1]))*wind[1]-sizeof(d);
		double dd[sizeof(d)+c];
		for(int i=0;i<sizeof(dd);i++) {
			if(i<sizeof(d)) { dd[i] = d[i]; } else { dd[i] = 0; }
		}
		
		double *m = mSum(dd,wind[x], fac);
		for(int y=0;y<size;y++) {
			mat[y][x] = m[y];
		}
	}
	
	int colc = sizeof(wind)/sizeof(int);
	for(int i=0;i<size;i++) {
		double tot =0;
		for(int j=0;j<colc;j++) {
			tot += mat[i][j];
		}
		vec[i] = tot/colc;
	}
	
	return(vec);
}

void cncp(double *dat, int *size) 
{
	double s=0;
	int ss = *size;
	double data[ss];
	for(int i=0;i<ss;i++) {
		s +=dat[0];
		data[i] = s;
		dat++;
	}
	double *d = diff(data, ss);
	for(int i=0;i<ss;i++) {
		dat--;
	}
	for(int i=0;i<ss;i++) {
		dat[0]=d[i];
		dat++;
	}

}




