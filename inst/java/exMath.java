import java.lang.Math.*;

public class exMath {

    public static double[] Absol(double[] values) {    
        double[] a = new double[values.length];
        for (int i =0; i<values.length ; i++)  a[i] = Math.abs(values[i]);    
        return a;    
    }
    
    public static float Round(float Rval, int Rpl) {
        float p = (float)Math.pow(10,Rpl);
        Rval = Rval * p; float tmp = Math.round(Rval);
        return (float)tmp/p;
    }

    public static double mean(double[] m) {
      double men = 0; int l = m.length;  
      for (int i = 0 ; i<m.length ; i++) men = men + m[i];
      return men / l;
   }
  
  	public static double mean(Object[] m) {
    	double men = 0; int l = m.length;    
      	for (int i = 0 ; i<m.length ; i++) men = men + Double.parseDouble(""+m[i]);
      	if (m.length==1) men =  Double.parseDouble(""+m[0]);  
    	return men / l;
  	}
  	
  	public static double median(double[] values) {    
      double[] copy = new double[values.length]; System.arraycopy(values, 0, copy, 0, copy.length);
      java.util.Arrays.sort(copy); int middle = copy.length/2;
      if (copy.length%2 == 1) return copy[middle];
	  else return (copy[middle-1] + copy[middle]) / 2.0;
	}
  	
  	public static double median(Object[] values1) {
  		double[] values = new double[values1.length];
  		for(int i=0;i<values.length;i++) values[i] = Double.parseDouble(""+values1[i]);
    	double[] copy = new double[values.length];System.arraycopy(values, 0, copy, 0, copy.length);
      	java.util.Arrays.sort(copy); int middle = copy.length/2;
      	if (copy.length%2 == 1) return copy[middle];
      	else if (copy.length==1) return(copy[0]);
      	else return (copy[middle-1] + copy[middle]) / 2.0;
	}

  public static double[] runMean(double[] m, int wind) {
      double[] retM = new double[m.length];
      for (int i=0;i<m.length-wind;i++) {
          double[] temp = new double[wind]; System.arraycopy(m, i, temp, 0, temp.length);
          double men = mean(temp);
          for (int j=i;j<i+wind;j++) retM[j] = men;
      }
      return(retM);
  }

  public static double[] runMedian(double[] m, int wind) {
      double[] retM = new double[m.length];
      for (int i=0;i<m.length-wind;i++) {
          double[] temp = new double[wind]; System.arraycopy(m, i, temp, 0, temp.length);
          double men = median(temp);
          for (int j=i;j<i+wind;j++) retM[j] = men;
      }
      return(retM);
  }

  public static double[] runMedianR(double[] m, int wind) {
      double[] retM = new double[m.length]; double[] rev = new double[m.length];
      int pin =0;
      for (int i=m.length-1;i>=0;i--) { rev[pin] = m[i]; pin++; }
      for (int i=0;i<rev.length-wind;i++) {
          double[] temp = new double[wind]; System.arraycopy(rev, i, temp, 0, temp.length);
          double men = median(temp);
          for (int j=i;j<i+wind;j++)  retM[j] = men;
      }
      pin =0;
      for (int i=m.length-1;i>=0;i--) { rev[pin] = retM[i]; pin++; }
      return(rev);
  }

   public static double Var(double[] data) {
        int n = data.length; int d = n-1;
        double[] xbar = new double[data.length]; 
        double m = mean(data); double var = 0;
        for (int i=0;i<n;i++) {
            double deviation = data[i] - m;
            var += Math.pow(deviation,2);
        }
        return var/d;
    }

   public static double quantile(double[] values, double quantile) {
        double[] copy = new double[values.length];System.arraycopy(values, 0, copy, 0, copy.length);
        java.util.Arrays.sort(copy); int index = (int) (copy.length * quantile);
        return copy[index];
    }

   public static double IQR(double[] values) {
       double q25 = quantile(values, 0.25); double q75 = quantile(values, 0.75);
       return q75 - q25;
   }

   public static double dLRs(double[] values) {
       double[] diff = new double[values.length -1];
       for (int i = 0; i<values.length -1;i++) {
           int k = i+1; diff[i] = values[i] - values[k];
       }
       double dLRs = IQR(diff) / 1.907745;  //  4 * qnorm((1+0.5) / 2) / sqrt(2)
       return dLRs;
   }

   public static double getPearsonCorrelation(double[] scores1,double[] scores2){
       double result = 0; double sum_sq_x = 0;
       double sum_sq_y = 0; double sum_coproduct = 0;
       double mean_x = scores1[0]; double mean_y = scores2[0];
        for(int i=2;i<scores1.length+1;i+=1){
            double sweep =Double.valueOf(i-1)/i; double delta_x = scores1[i-1]-mean_x; double delta_y = scores2[i-1]-mean_y;
            sum_sq_x += delta_x * delta_x * sweep; sum_sq_y += delta_y * delta_y * sweep;
            sum_coproduct += delta_x * delta_y * sweep;
            mean_x += delta_x / i; mean_y += delta_y / i;
        }
       double pop_sd_x = (double) Math.sqrt(sum_sq_x/scores1.length);
       double pop_sd_y = (double) Math.sqrt(sum_sq_y/scores1.length);
       double cov_x_y = sum_coproduct / scores1.length;
       result = cov_x_y / (pop_sd_x*pop_sd_y);
       return result;
   }


   public static double[] LinRes(double[] x, double[] y) {
       int MAXN = 10000; int n = 0;
       double[] res = new double[x.length]; double sumx = 0.0, sumy = 0.0, sumx2 = 0.0;
        if (x.length != y.length) System.exit(0);
        for (int i =0 ; i< x.length ; i++) {
            sumx  += x[n]; sumx2 += x[n] * x[n]; sumy  += y[n]; n++;
        }
        double xbar = sumx / n; double ybar = sumy / n;
        double xxbar = 0; double yybar = 0; double xybar = 0;
        for (int i = 0; i < n; i++) {
            xxbar += (x[i] - xbar) * (x[i] - xbar); yybar += (y[i] - ybar) * (y[i] - ybar); xybar += (x[i] - xbar) * (y[i] - ybar);
        }
        double beta1 = xybar / xxbar;  double beta0 = ybar - beta1 * xbar;
        int df = n - 2; double rss = 0; double ssr = 0;     
        for (int i = 0; i < n; i++) {
            double fit = beta1*x[i] + beta0;
            rss += (fit - y[i]) * (fit - y[i]); ssr += (fit - ybar) * (fit - ybar);
            res[i] = fit - y[i]; res[i] = -res[i];
        }
        return res;
    }

    public static double[] LinResSlope(double[] x, double[] y) {
        int MAXN = 10000; int n = 0;
        double[] res = new double[x.length]; double sumx = 0.0, sumy = 0.0, sumx2 = 0.0;
        if (x.length != y.length) System.exit(0);
        for (int i =0 ; i< x.length ; i++) {
            sumx  += x[n]; sumx2 += x[n] * x[n]; sumy  += y[n]; n++;
        }
        double xbar = sumx / n; double ybar = sumy / n;
        double xxbar = 0.0, xybar = 0.0;
        for (int i = 0; i < n; i++) {
            xxbar += (x[i] - xbar) * (x[i] - xbar); xybar += (x[i] - xbar) * (y[i] - ybar);
        }
        double beta1 = xybar / xxbar;
        for (int i = 0; i < n; i++) res[i] = x[i] / beta1;
        return res;
    }

	public static double deci(double x, double y, int p, int q) {		
	return( ((2*y)-(3*x)) *(q/(2*p)+0.5) );
	}
   

   	public double pscale(double[] x, double fac, int tnum, int minEst) {
		double dl = dLRs(x); double rp = quantile(x, 0.68);
		double L = ((Math.abs(rp-dl) + Math.abs(dl-rp)) / 2)*Math.sqrt(tnum/(minEst/tnum)*fac);
	return(fac*(fac*L+1));
   	}
	
}