#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>
#include        <iostream>
#include        "Library/TNT/tnt.h"
#include        "Library/MyFun/myFun.h"
using namespace std;
using namespace TNT; 

extern "C" {

static void    SubFun1(double &newInvV, double &newNu, double &newTmp,
                 double oldInvV, double oldNu, double ds, double obsds);
static void     BcmixForward(Matrix<double> &forNu,
                  Matrix<double> &forInvV, Matrix<double> &forQ,
                  Vector<double> &forP, Matrix<int> &forFilter,
                  Matrix<double> &forLog, Vector<double> obs,
                  double p, double a, double b, double mu, double v,
                  double sigma2, int K, int M);
static void     BcmixBackward(Matrix<double> &backNu,
                  Matrix<double> &backInvV, Matrix<double> &backQ,
                  Vector<double> &backP, Matrix<int> &backFilter,
                  Matrix<double> &backLog, Vector<double> obs,
                  double p, double a, double b, double mu, double v,
                  double sigma2, int K, int M);
static void   BcmixSmooth(Vector<double> &estPara, Vector<double> &estBS,
                Vector<double> obs, double p, double a,
                double b, double mu, double v, double sigma2, int K,
                int M);


/***************************************************************************
        This file contains all the code for mean shift models with
constant variance. The transition matrix
        (1-p    p/2     p/2
           c      a      b
           c      b      a)

    Note that to be consistent with the subscript of symbols in the 
paper, the index in the code starts from 1 (instead of 0 in the C/C++ 
convection). hence the input obsArray in the code becomes [0, obs].

***************************************************************************/


/* estBS = estimates of post prob of baseline state */
static void   BcmixSmooth(Vector<double> &estPara, Vector<double> &estBS,
        Vector<double> obs, double p, double a, 
        double b, double mu, double v, double sigma2, int K, 
        int M)
{  
  int       N = obs.size()-1;
  double        c=1-a-b, ds = 1/sigma2, iv = 1/v, muiv = mu/v,
                pstar, obsds, sumvec, tmp00=(mu*mu/v+log(v))/2;
  double        alpha, tmpa, tmpb=a/p*exp(tmp00),  total;
  Vector<double>   forP(N+1), backP(N+1);
  Matrix<int>      forFilter(N+1,K+1), backFilter(N+1,K+1);
  Matrix<double>   forNu(N+1,K+1), forInvV(N+1,K+1), forQ(N+1,K+1),
    forLog(N+1,K+1), backNu(N+1,K+1), backInvV(N+1,K+1), 
    backQ(N+1,K+1), backLog(N+1,K+1);

  BcmixForward(forNu, forInvV, forQ, forP, forFilter, forLog, obs,
    p, a, b, mu, v, sigma2, K, M);
  BcmixBackward(backNu, backInvV, backQ, backP, backFilter,
    backLog, obs, p, a, b, mu, v, sigma2, K, M);
  for (int t=1; t<N; t++) { estPara[t]=0;
    total = alpha = forP[t]*(c+(1-p-c)*backP[t+1])/c;
    tmpa = (b+(p-b)*backP[t+1])/p;
    for (int i=1;i<=min(t,K);i++)
      estPara[t] += forQ[t][i] * forNu[t][i]/forInvV[t][i];
    estPara[t] *= tmpa;     total += (1-forP[t])*tmpa;
    for (int i=1;i<=min(t,K);i++)  for (int j=1;j<=min(N+1-t,K);j++) {
      double    nuij, ivij, logij, incre;
      ivij = forInvV[t][i]+backInvV[t+1][j]-iv;
      nuij = forNu[t][i]+backNu[t+1][j]-muiv;
      logij= (nuij*nuij/ivij - log(ivij))*0.5;
      total += incre = tmpb *forQ[t][i]*backQ[t+1][j] * exp( logij
    -forLog[t][i] - backLog[t+1][j]);
      estPara[t] += incre * nuij/ivij;
    }
    estPara[t] /= total;    estBS[t] = alpha/total;
  }
  estPara[N]=0;         estBS[N] = forP[N];
  for (int i=1;i<=K;i++)
    estPara[N]+=forQ[N][i]*forNu[N][i]/forInvV[N][i];
}


static void     BcmixBackward(Matrix<double> &backNu,
                  Matrix<double> &backInvV, Matrix<double> &backQ,
                  Vector<double> &backP, Matrix<int> &backFilter,
                  Matrix<double> &backLog, Vector<double> obs,
                  double p, double a, double b, double mu, double v,
                  double sigma2, int K, int M)
{
  int           N = obs.size()-1, s;
  double        c=1-a-b, ds = 1/sigma2, iv = 1/v, muiv = mu/v,
                pstar, obsds, sumvec, tmp00=(mu*mu/v+log(v))/2;
  Vector<double>        backQstar(K+2);
  Vector<int>           tmpBackFil(K+2);
  Matrix<double>        tmpPost(K+2,4); //[][1]:InvV,[][2]:Nu,[][3]:Log
  
  backQ[N][1]=p/(p+c);           backP[N]=1-backQ[N][1];
  SubFun1(backInvV[N][1],backNu[N][1],backLog[N][1],iv,muiv,ds,obs[N]*ds);
  backFilter[N][1] = N;
  for (int t=N-1; t>=N-K+1; t--) {
    obsds=obs[t]*ds;    sumvec=pstar=(1-p-c)*backP[t+1]+c;    s=N+1-t;
    SubFun1(tmpPost[s][1],tmpPost[s][2],tmpPost[s][3],iv,muiv,ds,obsds);
    sumvec+=backQstar[s]=(b+(p-b)*backP[t+1])*exp(tmpPost[s][3]-tmp00);
    for (int i=1; i<s; i++) {
      SubFun1(tmpPost[i][1],tmpPost[i][2],tmpPost[i][3],
        backInvV[t+1][i], backNu[t+1][i], ds, obsds);
      sumvec += backQstar[i] = a*backQ[t+1][i]* exp(tmpPost[i][3]
    - backLog[t+1][i]);
    }
    backP[t] = pstar/sumvec;
    for (int i=1; i<=s; i++)  {
      backInvV[t][i] = tmpPost[i][1];    backNu[t][i] = tmpPost[i][2];
      backLog[t][i] = tmpPost[i][3];     backQ[t][i] = backQstar[i]/sumvec;
      backFilter[t][i] = N+1-i;
    }
  }

  for (int t=N-K; t>=1; t--) {
    obsds=obs[t]*ds;    sumvec=pstar=(1-p-c)*backP[t+1]+c;    s=N+1-t;
    SubFun1(tmpPost[K+1][1],tmpPost[K+1][2],tmpPost[K+1][3],iv,muiv,ds,obsds);
    sumvec+=backQstar[K+1]=(b+(p-b)*backP[t+1])*exp(tmpPost[K+1][3]-tmp00);
    tmpBackFil[K+1]=t;
    for (int i=1; i<K+1; i++) {
      SubFun1(tmpPost[i][1], tmpPost[i][2], tmpPost[i][3],
    backInvV[t+1][i], backNu[t+1][i], ds, obsds);
      sumvec += backQstar[i] = a*backQ[t+1][i]* exp(tmpPost[i][3]
    - backLog[t+1][i]);
      tmpBackFil[i] = backFilter[t+1][i];
    }
    backP[t] = pstar/sumvec;     int  minind = 1;
    for (int i=1;i<=K+1-M;i++) if (backQstar[minind]>backQstar[i]) minind=i;
    double delta=(1-backP[t])/(sumvec-pstar-backQstar[minind]);
    for (int i=1; i<minind; i++) {
      backInvV[t][i]=tmpPost[i][1];      backNu[t][i]=tmpPost[i][2];
      backLog[t][i]=tmpPost[i][3];       backFilter[t][i]=tmpBackFil[i];
      backQ[t][i]=backQstar[i]*delta;
    }
    for (int i=minind+1; i<=K+1; i++) {
      backInvV[t][i-1]=tmpPost[i][1];    backNu[t][i-1]=tmpPost[i][2];
      backLog[t][i-1]=tmpPost[i][3];     backFilter[t][i-1]=tmpBackFil[i];
      backQ[t][i-1]=backQstar[i]*delta;
    }
  }
}

static void BcmixForward(Matrix<double> &forNu, 
          Matrix<double> &forInvV, Matrix<double> &forQ,
          Vector<double> &forP, Matrix<int> &forFilter,
          Matrix<double> &forLog, Vector<double> obs, 
          double p, double a, double b, double mu, double v, 
          double sigma2, int K, int M)
{
  int           N = obs.size()-1;
  double        c=1-a-b, ds = 1/sigma2, iv = 1/v, muiv = mu/v,
                pstar, obsds, sumvec, tmp00=(mu*mu/v+log(v))/2;
  Vector<double>    forQstar(K+2);
  Vector<int>       tmpForFil(K+2);
  Matrix<double>    tmpPost(K+2,4); //[][1]:InvV,[][2]:Nu,[][3]:Log

  forQ[1][1]=p/(p+c);       forP[1]=1-forQ[1][1];
  SubFun1(forInvV[1][1],forNu[1][1],forLog[1][1],iv,muiv,ds,obs[1]*ds);
  forFilter[1][1] = 1;
  for (int t=2; t<=K; t++) {
    obsds=obs[t]*ds;    sumvec=pstar=(1-p-c)*forP[t-1]+c;
    SubFun1(tmpPost[t][1],tmpPost[t][2],tmpPost[t][3],iv,muiv,ds,obsds);
    sumvec+=forQstar[t]=(b+(p-b)*forP[t-1])*exp(tmpPost[t][3]-tmp00);
    for (int i=1; i<t; i++) {
      SubFun1(tmpPost[i][1],tmpPost[i][2],tmpPost[i][3],
    forInvV[t-1][i], forNu[t-1][i], ds, obsds);
      sumvec+=forQstar[i]=a*forQ[t-1][i]*exp(tmpPost[i][3]-forLog[t-1][i]);
    }
    forP[t] = pstar/sumvec;
    for (int i=1; i<=t; i++)  {
      forInvV[t][i] = tmpPost[i][1];    forNu[t][i] = tmpPost[i][2];
      forLog[t][i] = tmpPost[i][3]; forQ[t][i] = forQstar[i]/sumvec;
      forFilter[t][i] = i;
    }
  }

  for (int t=K+1; t<N+1; t++) {
    obsds=obs[t]*ds;    sumvec=pstar=(1-p-c)*forP[t-1]+c;
    SubFun1(tmpPost[K+1][1],tmpPost[K+1][2],tmpPost[K+1][3],iv,muiv,ds,obsds);
    sumvec+=forQstar[K+1]=(b+(p-b)*forP[t-1])*exp(tmpPost[K+1][3]-tmp00);
    tmpForFil[K+1]=t;
    for (int i=1; i<K+1; i++) {
      SubFun1(tmpPost[i][1],tmpPost[i][2],tmpPost[i][3],
        forInvV[t-1][i], forNu[t-1][i], ds, obsds);
      sumvec+=forQstar[i]=a*forQ[t-1][i]*exp(tmpPost[i][3]-forLog[t-1][i]);
      tmpForFil[i] = forFilter[t-1][i];
    }
    forP[t] = pstar/sumvec; int  minind = 1;
    for (int i=1;i<=K+1-M;i++) if (forQstar[minind]>forQstar[i]) minind=i;
    double delta=(1-forP[t])/(sumvec-pstar-forQstar[minind]);
    for (int i=1; i<minind; i++) {
      forInvV[t][i]=tmpPost[i][1];  forNu[t][i]=tmpPost[i][2];
      forLog[t][i]=tmpPost[i][3];   forFilter[t][i]=tmpForFil[i];
      forQ[t][i]=forQstar[i]*delta;
    }
    for (int i=minind+1; i<=K+1; i++) {
      forInvV[t][i-1]=tmpPost[i][1];    forNu[t][i-1]=tmpPost[i][2];
      forLog[t][i-1]=tmpPost[i][3]; forFilter[t][i-1]=tmpForFil[i];
      forQ[t][i-1]=forQstar[i]*delta;
    }
  }
}

static void    SubFun1(double &newInvV, double &newNu, double &newTmp,
         double oldInvV, double oldNu, double ds, double obsds)
{
  newInvV = oldInvV + ds;
  newNu = oldNu + obsds;
  newTmp = ( newNu*newNu/newInvV - log(newInvV) )*0.5;
}


void cppBcmixSmooth (double *obs, int *n, double *p, double *a, 
         double *b, double *mu, double *v, double *sigma2, 
         int *classify, int *K, int *M, double *sig, double *non0prob) {
     

    int N=*n;
    Vector<double>   obsvec(N+1), estPara(N+1), estBS(N+1);

    for (int i=N; i>0; i--) obsvec[i] = obs[i-1];  
    
    BcmixSmooth(estPara, estBS, obsvec, *p, *a, *b, *mu, *v, *sigma2, *K, *M);
    
    for (int i=0; i < N; i++) sig[i] = estPara[i+1]; 
    for (int i=0; i < N; i++) non0prob[i] = estBS[i+1];
}

void trymat(double *mat1, int *d1, int *d2, double *mat2){

    int D1 = *d1;
    int D2 = *d2;
    Matrix<double> mat1mat(D1,D2);
    int count=0;
    for(int i=0; i<D1; i++) 
        for(int j=0; j<D2; j++){
            mat1mat[i][j] = mat1[count];
            count++;
    }
    
    string filename = "mat1mat.txt";
    ofstream fout(filename.c_str());
    fout << mat1mat;
    fout.close();

    count=0;
    for(int i=0; i<D1; i++) 
        for(int j=0; j<D2; j++){
            mat2[count] = mat1mat[i][j];
            count++;
    }
}


}    /* extern "C" */

