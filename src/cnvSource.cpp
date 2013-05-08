#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>
using namespace std;
 
extern "C" {
       
class arrCGH {
      /**********************************************
          MAIN FUNCTIONS:
          
          void selectHyper(int);
          void getSmoothInfer();
          void getBcmixSmoothClass(int, int);
          void simOneSeq();
      ***********************************************/
      
      public:
      int n, K, M;
      double *obs, p, a, b, c, mu, v, sigma2;
      vector<double> expSig, expSig2, expSd, non0prob, 
                     expProbSig, expProbSig2, probCh, state0prob;
      vector<vector<double> > cmat, probseq;
      arrCGH(double *, int, double, double, double, double, 
         double, double, int, int, string);
           
      void selectHyper(int);
      void getSmoothInfer();
      void getBcmixSmoothClass(int, int);
      void simOneSeq();
      
      /**********************************************
          UTILITY FUNCTIONS:
          
          void getPostDist();
          void getForward();
          void getBackward();
          void getBcmixSmooth();
          void getBcmixForward();
          void getBcmixBackward();
          void smProbSeq(int, double);
          void emUpdateHyper();
      **********************************/
 
      private:
      double plogpi00, emTargetFun;
      vector<double> qFor, qBack, qBcmFor, qBcmBack, 
                     qitBcmStarFor, qitBcmStarBack;
      vector<vector<double> > pmu, pcov, psig2, plogpi, 
                              qitStarFor, qitStarBack, 
                              qitFor, qitBack,
                              betaij, betait, 
                              qitBcmFor, qitBcmBack;
      vector<vector<int> > qIndexBcmFor, qIndexBcmBack;
      void getPostDist();
      void getForward();
      void getBackward();
      void getBcmixSmooth();
      void getBcmixForward();
      void getBcmixBackward();
      void smProbSeq(int, double);
      void emUpdateHyper();      
};  /* class arrCGH */
 
arrCGH::arrCGH(double *obs1, int n1, double p1, double a1, double b1, 
        double mu1, double v1, double sigma21, int K1, int M1, string func) {
      
      obs = obs1; n = n1; p = p1; a = a1; b = b1; 
      c =  1 - a - b; mu = mu1; v = v1; sigma2 = sigma21;   
      
      K = K1; M = M1;  
      plogpi00 = 0; emTargetFun = 0;
           
      pmu.resize(n1); psig2.resize(n1); 
      pcov.resize(n1); plogpi.resize(n1); 
 
      if (func == "selHyp" || func == "getBcmSmCl") {
          
        qitBcmFor.resize(n1); qitBcmBack.resize(n1);
        qIndexBcmFor.resize(n1); qIndexBcmBack.resize(n1);
        probseq.resize(n1);
        
        betait.resize(n1);
 
        for (int i = 0; i < n1; i++) {
            pmu[i].resize(i+1, 0); psig2[i].resize(i+1, 0);
            pcov[i].resize(i+1, 0); plogpi[i].resize(i+1, 0); 
    
            qitBcmFor[i].resize(K, 0); qitBcmBack[i].resize(K, 0); 
            qIndexBcmFor[i].resize(K, -1); qIndexBcmBack[i].resize(K, -1); 
            probseq[i].resize(5, 0);
            
            betait[i].resize(K, 0);
          }
    
        qBcmFor.resize(n1, 0); qBcmBack.resize(n1, 0);
        qitBcmStarFor.resize(K+1, 0); qitBcmStarBack.resize(K+1, 0); 
      }     
 
      if (func == "getSmInfer") {
            
        qitFor.resize(n1); qitBack.resize(n1); cmat.resize(n1);
        betaij.resize(n1); betait.resize(n1); qitStarFor.resize(n1);
        qitStarBack.resize(n1); 
 
        for (int i = 0; i < n1; i++) {
         pmu[i].resize(i+1, 0); psig2[i].resize(i+1, 0);
         pcov[i].resize(i+1, 0); plogpi[i].resize(i+1, 0); 
    
         qitFor[i].resize(i+1, 0); qitBack[i].resize(i+1, 0); 
         cmat[i].resize(i+1, 0); betaij[i].resize(i+1, 0); 
         betait[i].resize(i+1, 0); qitStarFor[i].resize(i+1, 0); 
         qitStarBack[i].resize(i+1, 0); 
        }
 
        qFor.resize(n1, 0); qBack.resize(n1, 0);
      }     
                  
      expSig.resize(n1, 0); expSig2.resize(n1, 0); 
      expSd.resize(n1, 0); non0prob.resize(n1, 0); 
      expProbSig.resize(n1, 0); expProbSig2.resize(n1, 0); 
      probCh.resize(n1, 0); state0prob.resize(n1, 0); 
}
     
void arrCGH::selectHyper(int iterMax) {      
     getBcmixSmooth();
       
     int iter = 1;
     while (iter < iterMax-1) {                 
       emUpdateHyper();           
       getBcmixSmooth();
       iter++;            
     }
}
 
void arrCGH::getSmoothInfer() {
     getPostDist();
     getForward();
     getBackward();
     
     double total, alpha, expSigSum, expSig2Sum, expSig2SumN = 0;
     for (int k = 0; k < n-1; k++) {
         total = 0; expSigSum = 0; expSig2Sum = 0;
         alpha = (1-qFor[k])*(1-p+(c-1+p)*qBack[k+1])/c;
         total += alpha;
         for (int i = 0; i <= k; i++) {
             betait[k][i] = qitFor[k][i]*(p+(b-p)*qBack[k+1])/p;
             total += betait[k][i];
             expSigSum += betait[k][i]*pmu[k][i];
             expSig2Sum += betait[k][i]*psig2[k][i];
             for (int j = k+1; j < n; j++) {
                 betaij[k][i] = qitFor[k][i]*qitBack[j][k+1]* exp(plogpi[j][i] 
                              - plogpi[j][k+1] - plogpi[k][i] + plogpi00) *a/p;
                     
                 total += betaij[k][i];
                 expSigSum += betaij[k][i]*pmu[j][i];
                 expSig2Sum += betaij[k][i]*psig2[j][i];
             }
             cmat[k][i] = betait[k][i]/total;
             cmat[n-1][i] = qitFor[n-1][i];
         }
         non0prob[k] = 1 - alpha/total;
         expSig[k] = expSigSum/total;
         expSd[k] = sqrt(expSig2Sum/total - pow(expSig[k], 2));
         expSig[n-1] += pmu[n-1][k]*qitFor[n-1][k];
         expSig2SumN += psig2[n-1][k]*qitFor[n-1][k];
     }
     non0prob[n-1] = qFor[n-1];
     cmat[n-1][n-1] = qitFor[n-1][n-1];
     expSig[n-1] += pmu[n-1][n-1]*qitFor[n-1][n-1];
     expSig2SumN += psig2[n-1][n-1]*qitFor[n-1][n-1];
     expSd[n-1] = sqrt(expSig2SumN - pow(expSig[n-1], 2));
}
 
void arrCGH::getBcmixSmoothClass(int classify, int randNum) {
     getPostDist();
     getBcmixForward();
     getBcmixBackward();
     
     vector<double> beta, betait, mixobjW, mixobjMu, mixobjCov, temp1, temp2;
     mixobjW.reserve(420); mixobjMu.reserve(420); mixobjCov.reserve(420); 
 
     char buffer[15];
     string filename = "getBcmixSmoothClassOutput";
     sprintf(buffer, "%i", randNum);
     filename = filename.append(buffer) + ".dat";
     ofstream fout(filename.c_str());
     fout.precision(15);
     
     string filename2 = "BcmixCmat";
     filename2 = filename2.append(buffer) + ".dat";
     ofstream fout2(filename2.c_str());
     fout2.precision(15);
     
     for (int k = 0; k < n-1; k++) {
         double betaSum = 0, betaValue, sum = 0, sumMu = 0;   
         int bValue, bCount = -1, fValue, fCount = -1;
         
         beta.clear(); betait.clear();
         mixobjW.clear(); mixobjMu.clear(); mixobjCov.clear(); 
         temp1.clear(); temp2.clear();
         for (int i = 0; i < K; i++) {
             fValue = qIndexBcmFor[k][i];
 
             if (fValue >= 0) { 
                fCount++;               
                
                betait.push_back(qitBcmFor[k][fCount]*(p+(b-p)*qBcmBack[k+1])/p);
                mixobjW.push_back(qitBcmFor[k][fCount]*(p+(b-p)*qBcmBack[k+1])/p);
                mixobjMu.push_back(pmu[k][fValue]);
                mixobjCov.push_back(pcov[k][fValue]);
 
                sum += mixobjMu[fCount]*mixobjW[fCount];
                sumMu += mixobjMu[fCount];
               
                for (int j = 0; j < K; j++) {
                    bValue = qIndexBcmBack[n-k-2][j];
                    if (bValue >= 0) { 
                       bCount++;              
                       betaValue = qitBcmBack[n-k-2][bCount]*qitBcmFor[k][fCount]*exp(plogpi[bValue][fValue] 
                        - plogpi[bValue][k+1] - plogpi[k][fValue] + plogpi00)*a/p;  
                       beta.push_back(betaValue);
                       betaSum += betaValue;
                                              
                       temp1.push_back(pmu[bValue][fValue]);
                       temp2.push_back(pcov[bValue][fValue]);
 
                   sum += betaValue*pmu[bValue][fValue];
                   sumMu += pmu[bValue][fValue];
                   }
                }  
                bCount = -1;
             }
         }
 
         mixobjW.insert(mixobjW.end(), beta.begin(), beta.end());
         mixobjMu.insert(mixobjMu.end(), temp1.begin(), temp1.end());
         mixobjCov.insert(mixobjCov.end(), temp2.begin(), temp2.end());
                     
         smProbSeq(k, betaSum);
         double pTotal = probseq[k][0] + probseq[k][1] + probseq[k][2] 
                   + probseq[k][3] + probseq[k][4],
             alpha = probseq[k][0] + probseq[k][1];
     
         state0prob[k] = alpha/pTotal;
         expSig[k] = sum/pTotal;
         
         for (unsigned int m = 0; m < betait.size(); m++) {
             fout2 << betait[m]/pTotal << "@";
         }
         fout2 << endl;
 
     if (classify > 0) {
         for (unsigned int i=0; i < mixobjMu.size()-1; i++) fout << mixobjMu[i] << ",";
         fout << mixobjMu[mixobjMu.size()-1] << "@";
         for (unsigned int i=0; i < mixobjCov.size()-1; i++) fout << mixobjCov[i] << ",";
         fout << mixobjCov[mixobjCov.size()-1] << "@";
         for (unsigned int i=0; i < mixobjW.size()-1; i++) fout << mixobjW[i] << ",";
         fout << mixobjW[mixobjW.size()-1] << endl; }
     }  
 
     expSig[n-1] = 0; 
     for (int i=0; i < K; i++)
       expSig[n-1] += pmu[n-1][qIndexBcmFor[n-1][i]]*qitBcmFor[n-1][i];
     
     fout.close();
     fout2.close();
}
 
void arrCGH::getPostDist() {
      for (int i = 0; i < n; i++) {
          double cumsum = 0;
          for (int j = i; j < n; j++) {
              cumsum += obs[j]; 
              pcov[j][i] = pow(1/v + (j-i+1)/sigma2, -1);
              pmu[j][i] = (mu/v + cumsum/sigma2)*pcov[j][i];
              plogpi[j][i] = 0.5*(pow(pmu[j][i], 2)/pcov[j][i] + log(pcov[j][i]));
              psig2[j][i] = pow(pmu[j][i], 2) + pcov[j][i];
              
              
          }
      }
      plogpi00 = 0.5*(pow(mu, 2)/v + log(v));
}
 
void arrCGH::getForward() { 
     //this function computes the Bayes forward filter 
     qFor[0] = p/(p + c);
     qitFor[0][0] = p/(p + c);
     
     for (int k = 1; k < n; k++) {
         double cumsum = 0;
         for (int j = 0; j < k; j++) {
            qitStarFor[k][j] = exp(log(a*qitFor[k-1][j]) +
                                 plogpi[k][j] - plogpi[k-1][j]);
            cumsum += qitStarFor[k][j];
         }
         
         qitStarFor[k][k] = exp(log(p+(b-p)*qFor[k-1]) 
                               + plogpi[k][k] - plogpi00);
         cumsum += qitStarFor[k][k];
         
         double total = (1-p+(c-1+p)*qFor[k-1]) + cumsum;
         qFor[k] = cumsum/total;
         
         for (int j = 0; j <= k; j++) {
             qitFor[k][j] = qitStarFor[k][j]/total;
         }
     }
}
 
void arrCGH::getBackward() {     
     //this function computes the Bayes backward filter
     qBack[n-1] = p/(p + c);
     qitBack[n-1][n-1] = p/(p + c);
     
     for (int k=n-2; k>=0; k--) {
         double cumsum = 0;
         qitStarBack[k][k] = exp(log(p+(b-p)*qBack[k+1]) 
                               + plogpi[k][k] - plogpi00);
         cumsum += qitStarBack[k][k];
         for (int j = k+1; j < n; j++) {
            qitStarBack[j][k] = exp(log(a*qitBack[j][k+1]) +
                                 plogpi[j][k] - plogpi[j][k+1]);
            cumsum += qitStarBack[j][k];
         }
         
         double total = (1-p+(c-1+p)*qBack[k+1]) + cumsum;
         qBack[k] = cumsum/total;
         
         for (int j = k; j < n; j++) {
             qitBack[j][k] = qitStarBack[j][k]/total;
         }
     }
}
 
void arrCGH::getBcmixSmooth() {
     getPostDist();
     getBcmixForward(); 
     getBcmixBackward();
     
     vector<double> logVec(5);
     logVec[0] = log(1-p); logVec[1] = log(p); logVec[2] = log(c);
     logVec[3] = log(b); logVec[4] = log(a);
     double pi = 2*acos(0);
     
     emTargetFun = 0;
     for (int k = 0; k < n-1; k++) {
         double beta = 0, betaSum = 0, tmp; 
         vector<double> sum(6, 0);     
         int bValue, bCount = -1, fValue, fCount = -1;
         
         for (int i = 0; i < K; i++) {
             fValue = qIndexBcmFor[k][i];
             if (fValue >= 0) { 
                fCount++; 
 
                sum[0] += pmu[k][fValue]*qitBcmFor[k][fCount];
                sum[1] += psig2[k][fValue]*qitBcmFor[k][fCount];
                   
                for (int j = 0; j < K; j++) {
                    bValue = qIndexBcmBack[n-k-2][j];
                    if (bValue >= 0) { 
                       bCount++;              
                       beta = qitBcmBack[n-k-2][bCount]*qitBcmFor[k][fCount]*exp(plogpi[bValue][fValue] 
                              - plogpi[bValue][k+1] - plogpi[k][fValue] + plogpi00)*a/p;  
                       betaSum += beta;
                        
                       sum[2] += beta*psig2[bValue][fValue];     
                       
                       if (fCount == 0) {
                          tmp = (1+(b/p-1)*qBcmFor[k])*qitBcmBack[n-k-2][bCount];
                          sum[3] += tmp*pmu[bValue][k+1];
                          sum[4] += tmp*(psig2[bValue][k+1] - 2*mu*pmu[bValue][k+1] + pow(mu, 2));
                          sum[5] += tmp;
                       }
                   }
                }  
                bCount = -1;
             }
         }
 
         smProbSeq(k, betaSum);
         double pTotal = probseq[k][0] + probseq[k][1] + probseq[k][2] 
                       + probseq[k][3] + probseq[k][4];      
         for (int i = 0; i < 5; i++) probseq[k][i] = probseq[k][i]/pTotal;
     
         non0prob[k] = 1 - probseq[k][0] - probseq[k][1];
         
         expSig[k] = (sum[0]*(p+(b-p)*qBcmBack[k+1])/p)/pTotal;
         expSig2[k] = (sum[1]*(p+(b-p)*qBcmBack[k+1])/p + sum[2])/pTotal; 
         expProbSig[k] = sum[3]/pTotal;
         expProbSig2[k] = sum[4]/pTotal;
         probCh[k] = sum[5]/pTotal;
                  
         emTargetFun += -(pow(obs[k],2) - 2*obs[k]*expSig[k] 
              + expSig2[k])/2/sigma2 + (probseq[k][0] + probseq[k][1] 
              + probseq[k][2] + probseq[k][3] + probseq[k][4])*logVec[k%5]
              - 0.5*expProbSig2[k]/v - 0.5*log(v)*probCh[k] 
              - 0.5*log(2*pi)*(probseq[k][1] + probseq[k][3]);
     }  
     
     expSig[n-1] = 0; expSig2[n-1] = 0;
     for (int i = 0; i < K; i++) {
         expSig[n-1] += pmu[n-1][qIndexBcmFor[n-1][i]]*qitBcmFor[n-1][i];
         expSig2[n-1] += psig2[n-1][qIndexBcmFor[n-1][i]]*qitBcmFor[n-1][i];
     }
     
     emTargetFun += -(pow(obs[n-1],2) - 2*obs[n-1]*expSig[n-1] 
         + expSig2[n-1])/2/sigma2 - n*log(sigma2)/2 + (probseq[n-1][0] 
         + probseq[n-1][1] + probseq[n-1][2] + probseq[n-1][3] 
         + probseq[n-1][4])*logVec[(n-1)%5] - 0.5*expProbSig2[n-1]/v 
         - 0.5*log(v)*probCh[n-1] - 0.5*log(2*pi)*(probseq[n-1][1] 
         + probseq[n-1][3]); 
}
 
void arrCGH::getBcmixForward() {
     //approximations of forward filters by BCMIX
     qBcmFor[0] = p/(p + c);
     qitBcmFor[0][0] = p/(p + c);
     qIndexBcmFor[0][0] = 0;
     
     for (int k = 1; k < K; k++) {
         double cumsum = 0;
         for (int j = 0; j < k; j++) {
            qitBcmStarFor[j] = exp(log(a*qitBcmFor[k-1][j]) +
                                 plogpi[k][j] - plogpi[k-1][j]);
            cumsum += qitBcmStarFor[j];
         }
         
         qitBcmStarFor[k] = exp(log(p+(b-p)*qBcmFor[k-1]) 
                               + plogpi[k][k] - plogpi00);
         cumsum += qitBcmStarFor[k];
 
         double total = (1-p+(c-1+p)*qBcmFor[k-1]) + cumsum;
         qBcmFor[k] = cumsum/total;
         
         for (int j = 0; j <= k; j++) {
             qitBcmFor[k][j] = qitBcmStarFor[j]/total;
             qIndexBcmFor[k][j] = j;
         }
     }
     
     for (int k = K ; k < n; k++) {
         double cumsum = 0;
         for (int j = 0; j < K; j++) {
            int index = qIndexBcmFor[k-1][j];
            qitBcmStarFor[j] = exp(log(a*qitBcmFor[k-1][j]) +
                                 plogpi[k][index] - plogpi[k-1][index]);
            cumsum += qitBcmStarFor[j];
         }
      
         qitBcmStarFor[K] = exp(log(p+(b-p)*qBcmFor[k-1]) 
                               + plogpi[k][k] - plogpi00);
         cumsum += qitBcmStarFor[K];
               
         double total = (1-p+(c-1+p)*qBcmFor[k-1]) + cumsum;
         qBcmFor[k] = cumsum/total;
         
         int minIndex = 0;
         for (int j = 0; j < M; j++) {
             if (qitBcmStarFor[j] < qitBcmStarFor[minIndex])
                minIndex = j;
         }
                
         int index = 0;
         for (int j = 0; j <= K; j++) {
             if (j != minIndex) {
                  qitBcmFor[k][index] 
                    = qitBcmStarFor[j]/total*cumsum/(cumsum - qitBcmStarFor[minIndex]);
                  qIndexBcmFor[k][index] = qIndexBcmFor[k-1][j];
                  index++;
             }
         }
         qIndexBcmFor[k][K-1] = k;
     }
}
 
void arrCGH::getBcmixBackward() {
     //approximations of backward filters by BCMIX
     qBcmBack[0] = p/(p + c);
     qitBcmBack[0][0] = p/(p + c);
     qIndexBcmBack[0][0] = n-1;
     
     for (int k = 1; k < K; k++) {
         double cumsum = 0;
         for (int j = 0; j < k; j++) {
            qitBcmStarBack[j] = exp(log(a*qitBcmBack[k-1][j]) +
                                 plogpi[n-j-1][n-k-1] - plogpi[n-j-1][n-k]);
            cumsum += qitBcmStarBack[j];
         }
    
         qitBcmStarBack[k] = exp(log(p+(b-p)*qBcmBack[k-1]) + 
                                    plogpi[n-k-1][n-k-1] - plogpi00);
         cumsum += qitBcmStarBack[k];
 
         double total = (1-p+(c-1+p)*qBcmBack[k-1]) + cumsum;
         qBcmBack[k] = cumsum/total;
         
         for (int j = 0; j <= k; j++) {
             qitBcmBack[k][j] = qitBcmStarBack[j]/total;
             qIndexBcmBack[k][j] = n-j-1;
         }    
     }
     
     for (int k = K ; k < n; k++) {
         double cumsum = 0;
         for (int j = 0; j < K; j++) {
            int index = qIndexBcmBack[k-1][j];
            qitBcmStarBack[j] = exp(log(a*qitBcmBack[k-1][j]) +
                                 plogpi[index][n-k-1] - plogpi[index][n-k]);
            cumsum += qitBcmStarBack[j];
         }
         
         qitBcmStarBack[K] = exp(log(p+(b-p)*qBcmBack[k-1]) 
                               + plogpi[n-k-1][n-k-1] - plogpi00);
         cumsum += qitBcmStarBack[K];
               
         double total = (1-p+(c-1+p)*qBcmBack[k-1]) + cumsum;
         qBcmBack[k] = cumsum/total;
         
         int minIndex = 0;
         for (int j = 0; j < M; j++) {
             if (qitBcmStarBack[j] < qitBcmStarBack[minIndex])
                minIndex = j;
         }
                
         int index = 0;
         for (int j = 0; j <= K; j++) {
             if (j != minIndex) {
                  qitBcmBack[k][index] 
                    = qitBcmStarBack[j]/total*cumsum/(cumsum - qitBcmStarBack[minIndex]);
                  qIndexBcmBack[k][index] = qIndexBcmBack[k-1][j];
                  index++;
             }
         }
         qIndexBcmBack[k][K-1] = n-k-1;
     }
}
 
void arrCGH::smProbSeq(int tt, double betaTot) {
     probseq[tt][0] = (1-qBcmFor[tt])*(1-p)*(1-qBcmBack[tt+1])/c;
     probseq[tt][1] = (1-qBcmFor[tt])*qBcmBack[tt+1];
     probseq[tt][2] = qBcmFor[tt]*(1-qBcmBack[tt+1]);
     probseq[tt][3] = qBcmFor[tt]*qBcmBack[tt+1]*b/p;
     probseq[tt][4] = betaTot;  
}
 
void arrCGH::emUpdateHyper() {
    //this functions computes the updated hyperparameters in EM procedure.
     vector<double> sum(9, 0);
     for (int i = 0; i < n; i++) {
         sum[0] += probseq[i][1];
         sum[1] += probseq[i][2];
         sum[2] += probseq[i][3];
         sum[3] += probseq[i][0] + probseq[i][1];
         sum[4] += probseq[i][2] + probseq[i][3] + probseq[i][4];
         sum[5] += pow(obs[i], 2) - 2*obs[i]*expSig[i] + expSig2[i];
         sum[6] += expProbSig[i];
         sum[7] += expProbSig2[i];
         sum[8] += probCh[i];
     }
                     
     p = sum[0]/sum[3]; b = sum[2]/sum[4]; c = sum[1]/sum[4];
     sigma2 = sum[5]/n; v = sum[7]/sum[8]; mu = sum[6]/sum[8];
     a = 1 - b - c;
}
 
void selectHyper (double *obs, int *n1, double *p0, 
      double *q0, int *K1, int *M1, int *iterMax, double *output) {
     
     int n, K, M;
     double p, q, a, b, mu, v, sigma2, sumX = 0, sumX2 = 0, 
            w, wsumX = 0, wsumX2 = 0;
     n = (*n1);
     p = (*p0);
     q = (*q0);
     a = 1 - q;
     b = q/2;
     K = *K1;
     M = *M1;
     
     int winsize = 50, len = n-winsize+1;;
     for (int i = 0; i < len; i++) {
         w = 0;
         for (int j = i; j <= (i+winsize-1); j++) {
             w += obs[j];
         }
         wsumX += w/winsize; 
         wsumX2 += pow(w/winsize, 2);
     }
     
     v = wsumX2/(len-1) - pow(wsumX, 2)/(len*(len-1));
     
     for (int i=0; i < n; i++) {
         sumX += obs[i];
         sumX2 += pow(obs[i], 2);
     }
     mu = sumX/n;
     
     sigma2 = sumX2/(n-1) - pow(sumX, 2)/(n*(n-1)) - v;
 
     arrCGH aCGH(obs, n, p, a, b, mu, v, sigma2, K, M, "selHyp");
     aCGH.selectHyper(*iterMax);
 
     output[0] = aCGH.p;
     output[1] = aCGH.a;
     output[2] = aCGH.b;
     output[3] = aCGH.mu;
     output[4] = aCGH.v;
     output[5] = aCGH.sigma2;
     
     /*ofstream fout("selectHyperOutput.dat");
     for (int i=0; i < n; i++) {
         fout << aCGH.expSig[i] << "," << aCGH.non0prob[i] << endl;
     }
     fout.close();*/
}
 
void getSmoothInfer (double *obs, int *n, double *p, double *a, double *b, 
          double *mu, double *v, double *sigma2, int *randNum) {
     
     int K = -1, M = -1;
     arrCGH aCGH(obs, *n, *p, *a, *b, *mu, *v, *sigma2, K, M, "getSmInfer");
     aCGH.getSmoothInfer();
     
     char buffer[15];
     string filename = "getSmoothInferOutput";
     sprintf(buffer, "%i", *randNum);
     filename = filename.append(buffer) + ".dat";
     ofstream fout(filename.c_str());
     fout.precision(15);
 
     for (int i=0; i < *n; i++) {
         fout << aCGH.expSig[i] << "," << aCGH.expSd[i] << ","
              << aCGH.non0prob[i];
         for (int j = 0; j < *n; j++) {
             if (j <= i) fout << "," << aCGH.cmat[i][j];
             else fout << ",";
         }
         fout << endl;
     }
     fout.close();
}
 
void getBcmixSmoothClass (double *obs, int *n, double *p, double *a, 
         double *b, double *mu, double *v, double *sigma2, 
         int *classify, int *K, int *M, int *randNum, double *output) {
     
     arrCGH aCGH(obs, *n, *p, *a, *b, *mu, *v, *sigma2, *K, *M, "getBcmSmCl");
     aCGH.getBcmixSmoothClass(*classify, *randNum);
     
     for (int i=0; i < *n; i++) 
       output[i] = aCGH.expSig[i];
}
} /* extern "C" */
