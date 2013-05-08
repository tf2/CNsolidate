/*
Author:	Tomas William Fitzgerald
All rights reserved.
*/

import java.io.*;
import java.util.*;
import java.lang.Math.*;
import java.util.ArrayList;


public class Wave {

	public static double factor = 0.5;
	public static double sized = 50;
	
	final static double[] g0a = {-0.0046, -0.0054, 0.0170, 0.0238, -0.1067, 0.0119, 0.5688, 0.7561, 0.2753, -0.1172, -0.0389, 0.0347, -0.0039, 0.0033};
    final static double[] g0b = {0.0033, -0.0039, 0.0347, -0.0389, -0.1172, 0.2753, 0.7561, 0.5688, 0.0119, -0.1067, 0.0238, 0.0170, -0.0054, -0.0046};
    final static double[] g1a = {-0.0033, -0.0039, -0.0347, -0.0389, 0.1172, 0.2753, -0.7561, 0.5688, -0.0119, -0.1067, -0.0238, 0.0170, 0.0054, -0.0046};
    final static double[] g1b = {-0.0046, 0.0054, 0.0170, -0.0238, -0.1067, -0.0119, 0.5688, -0.7561, 0.2753, 0.1172, -0.0389, -0.0347, -0.0039, -0.0033};	
    final static double[] g0o = {0.0001, 0, -0.0013, -0.0019, 0.0072, 0.0239, -0.0556, -0.0517, 0.2998, 0.5594, 0.2998, -0.0517, -0.0556, 0.0239, 0.0072, -0.0019, -0.0013, 0, 0.0001};
    final static double[] g1o = {-0.0018, 0, 0.0223, 0.0469, -0.0482, -0.2969, 0.5555, -0.2969, -0.0482, 0.0469, 0.0223, 0, -0.0018};
		
    final static double[] h0a = {0.0033, -0.0039, 0.0347, -0.0389, -0.1172, 0.2753, 0.7561, 0.5688, 0.0119, -0.1067, 0.0238, 0.0170, -0.0054, -0.0046};	
    final static double[] h0b = {-0.0046, -0.0054, 0.0170, 0.0238, -0.1067, 0.0119, 0.5688, 0.7561, 0.2753, -0.1172, -0.0389, 0.0347, -0.0039, 0.0033};	
    final static double[] h1a = {-0.0046, 0.0054, 0.0170, -0.0238, -0.1067, -0.0119, 0.5688, -0.7561, 0.2753, 0.1172, -0.0389, -0.0347, -0.0039, -0.0033};
    final static double[] h1b = {-0.0033, -0.0039, -0.0347, -0.0389, 0.1172, 0.2753, -0.7561, 0.5688, -0.0119, -0.1067, -0.0238, 0.0170, 0.0054, -0.0046};
	final static double[] h0o = {-0.0018, 0, 0.0223, -0.0469, -0.0482, 0.2969, 0.5555, 0.2969, -0.0482, -0.0469, 0.0223, 0, -0.0018};
    final static double[] h1o = {-0.0001, 0, 0.0013, -0.0019, -0.0072, 0.0239, 0.0556, -0.0517, -0.2998, 0.5594, -0.2998, -0.0517, 0.0556, 0.0239, -0.0072, -0.0019, 0.0013, 0, -0.0001};
	
	public Wave() {
	}
	
	public static ArrayList[] Input(String filename) {		
     	ArrayList[] data = new ArrayList[4];
		for (int i=0;i<data.length;i++) data[i] = new ArrayList();
		filename.trim(); String inLine; BufferedReader infile = null;	
		try {
			infile = new BufferedReader(new FileReader (filename));
			while ((inLine=infile.readLine()) != null) {	
				String[] Ldata = inLine.split("\t");
				data[0].add(Ldata[0]); data[1].add(Ldata[1]); data[2].add(Ldata[2]); data[3].add(Ldata[3]);				
			}
		} catch (FileNotFoundException ex) { System.out.println("File not found: " + filename); }
		  catch (IOException ex) { System.out.println(ex.getMessage()); }
		  finally { try { if (infile != null) infile.close(); } catch (IOException ex) { System.out.println(ex.getMessage()); }
		}	
	return data;	
    }
	
	public static double[] Absol(double[] values) {      
        double[] a = new double[values.length];	
        for (int i =0; i<values.length ; i++)  a[i] = Math.abs(values[i]); 
    return a;    
    }
	
	public static double mean(Object[] nn) {	
		double[] m = new double[nn.length];
		for(int i=0;i<m.length;i++) m[i]=Double.parseDouble(""+nn[i]);	
		double men = 0; int l = m.length;	
		for (int i = 0 ; i<m.length ; i++) men = men + m[i];     	
	return men / l;
	}
	
	public static double median(double[] values) { 
    	double[] copy = new double[values.length];
      	System.arraycopy(values, 0, copy, 0, copy.length); java.util.Arrays.sort(copy);
      	int middle = copy.length/2;
      	if (copy.length%2 == 1) return copy[middle]; else return (copy[middle-1] + copy[middle]) / 2.0;
	}
	
	public static double[] runMedian(double[] m, int wind) {
      double[] retM = new double[m.length];
      	for (int i=0;i<m.length-wind;i++) {
        	double[] temp = new double[wind];
          	System.arraycopy(m, i, temp, 0, temp.length);
          	double men = median(temp);
          	for (int j=i;j<i+wind;j++) {
              	retM[j] = men;
          	}
      	}
    	return(retM);
  	}
	
	public static double quantile(double[] values, double quantile) {
    	double[] copy = new double[values.length];
        System.arraycopy(values, 0, copy, 0, copy.length); java.util.Arrays.sort(copy);
        int index = (int) (copy.length * quantile);
    return copy[index];
    }
       
	public static double IQR(double[] values) {
   		double q25 = quantile(values, 0.25); double q75 = quantile(values, 0.75);
   	return q75 - q25;
   	}

   	public static double[] diff(double[] values) {
   		double[] diff = new double[values.length -1];
       	for (int i = 0; i<values.length -1;i++) {
           int k = i+1;diff[i] = values[i] - values[k];
       	}
   	return diff;
   	}

   	public static double dLRs(double[] values) {
       	double[] diff =diff(values); double dLRs = IQR(diff) / 1.907745;
   	return dLRs;
   	}
   
   public static double[] rtab(double[] p1, double[] p2) {  
		double[] pp = new double[p1.length]; double ar = quantile(Absol(p2),0.68); double r = quantile(Absol(p1),0.68)*factor;
		int size = (int)sized;
		for(int i=1;i<p2.length-1;i++) {
			double ptmd = 0;
			double[] temp = new double[5]; int pin = i; if(pin>=p2.length-4) { pin = pin-3; }
          	System.arraycopy(p2, pin, temp, 0, temp.length); double tmd = median(temp);
          	temp = new double[25]; int pin2 = i; if(pin2>=p2.length-(24)) { pin2 = pin2-23; }
          	System.arraycopy(p2, pin2, temp, 0, temp.length); double tmd11 = median(temp);
          	temp = new double[size]; pin2 = i; if(pin2>=p2.length-(size-1)) { pin2 = pin2-(size-2); }
          	System.arraycopy(p2, pin2, temp, 0, temp.length); double tmd1 = median(temp);
          	temp = new double[size]; pin2 = i-size; if(pin2<0) { pin2 =0; }
          	System.arraycopy(p2, pin2, temp, 0, temp.length); double rtmd1 = median(temp);
          	temp = new double[5]; pin2 = i; if(pin2>=p2.length-4) { pin2 = pin2-3; }
          	System.arraycopy(p2, pin2, temp, 0, temp.length); double tmd15 = median(temp);
          	temp = new double[5]; pin2 = i-5; if(pin2<0) { pin2 =0; }
          	System.arraycopy(p2, pin2, temp, 0, temp.length); double rtmd15 = median(temp);
          	for(int j=0;j<temp.length;j++) { 
          		temp[j] = temp[j]-tmd1; 
          		if(Math.abs(temp[j])>ar) { temp[j] = tmd; }
          	}
          	double rr = quantile(Absol(temp),0.68); System.arraycopy(p1, pin2, temp, 0, temp.length);
			double tmd2 = median(temp); double las = 0;	if(Math.abs(tmd2) < ar) { ptmd=tmd2; }
			if(Math.abs(p2[i])<r) {
 				pp[i] = p2[i]-p1[i];	
 				if(Math.abs(p1[i]-tmd2) < Math.abs(rr*Math.abs(tmd11*(factor)))) { las = p1[i]; }
 				//if(Math.abs(p1[i]-tmd2) < 2*Math.abs(rr*Math.abs(tmd11*(factor)))) { las = p1[i]; }
				if(Math.abs(p1[i]-tmd2) > Math.abs(rr*Math.abs(tmd11*factor))) { pp[i] = p2[i]-las; }		
				if(Math.abs(tmd1)>r & Math.abs(tmd2)>rr) { pp[i] = p2[i]-las; }
				if(Math.abs(p2[i])>rr) { pp[i] = p2[i]-las; }	
				if(Math.abs(p2[i-1])>rr | Math.abs(p2[i+1])>rr) { pp[i] = p2[i]-las; }
			} else {	
				//las = p1[i];
				//if(Math.abs(p1[i]) < Math.abs(rr*Math.abs(tmd1*(factor)))) { las = p1[i]; }
				if(Math.abs(p1[i]) < Math.abs(tmd) & Math.abs(tmd) < Math.abs(ar) ) { las = p1[i]; }
				pp[i] = p2[i]-las; 
			}
			// Three hacks... the kind that just work... remove at your own risk!
			if(Math.abs(tmd1)>(ar*(factor/2)) | Math.abs(rtmd1)>(ar*(factor/2)) | Math.abs(tmd15)>(ar*factor) | Math.abs(rtmd15)>(ar*factor)) { if(Math.abs(tmd11)<ar*factor) { pp[i]=(p2[i]-p1[i])+tmd11; } else { pp[i]=p2[i]-ptmd; } }
			if(Math.abs(pp[i]-p2[i])>rr) { if(Math.abs(las)<ar*factor) { pp[i] = p2[i]-las; } else { pp[i] = p2[i]-ptmd; } }
			if(pp[i]>0 & p2[i]<0 | pp[i]<0 & p2[i]>0) { if(Math.abs(pp[i])>ar*3) { pp[i]=p2[i]-ptmd; } }	
		}
   	return(pp);
   }
 	    
    public static double[] extend(double[] data) {	
        double[] newdata = new double[data.length+2]; newdata[0] = data[0];
        for (int i=0;i<data.length;i++) newdata[i+1] = data[i];
        newdata[newdata.length-1] = data[data.length-1];
    return newdata;
    }
    
    public static double[] convolve(double[] data, double[] operator, double[] output){	
        int dataLen = data.length; int operatorLen = operator.length;	
        for(int i = 0;i < dataLen-operatorLen+1;i++){
            output[i] = 0; for(int j = operatorLen-1;j >= 0;j--) output[i] += data[i+j]*operator[j];
        }
    return output;
    }
   	
   	public static double[] Reflect(double[] dataX, double minx, double maxx) {	
        double[] dataY = new double[dataX.length];  
        for (int i=0;i<dataX.length;i++) {		
            if (dataX[i]>maxx) dataY[i] = (2* maxx) - dataX[i];
            else dataY[i] = dataX[i];
        }
        for (int i=0;i<dataY.length;i++) {		
            if (dataY[i]<minx) dataY[i] = (2* minx) - dataY[i];
            else dataY[i] = dataX[i];
            if (dataX[i]>maxx) dataY[i] = (2* maxx) - dataY[i];
            else dataY[i] = dataY[i];	
        }
    return dataY;
    }
   	
   	public static double[] dtWaveIfm(double[] Yl, double[] Yh, int nLevels) {		
		double[] Lo = Yl;
		for (int i =0;i<nLevels;i++) Lo = coliFilt(Yl, g0b, g0a);	
    return Lo;
    }
	
    public static double[] dtWaveXfm(double[] data, int nLevels) {		
        int L = data.length; boolean odd = false;
        if (L % 2 == 1 ) { odd=true; L++; }	
		double[] data2 = new double[L];
		for(int i=0;i<L-1;i++)  data2[i] = data[i];
		if(odd) data2[L-1] = 0;		
        double[] Hi = colFilter(data2, h1o); double[] Lo = colFilter(data2, h0o);    
        for (int i=0;i<nLevels;i++) {	
            if (Lo.length<10) break;	
            if (Lo.length % 4 != 0 ) Lo = extend(Lo);  	
            Hi = coldFilt(Lo,h1b,h1a); Lo = coldFilt(Lo,h0b,h0a);
			Lo = dtWaveIfm(Lo, Hi, nLevels);		
        }      
        for(int i=Lo.length-1;i<0;i++) {
        	Lo[i]=Lo[i-1];
        }
    return Lo;	
    }
		
    public static double[] colFilter(double[] data, double[] filter) {	
        int r = data.length; int m = filter.length; int m2 = (int) Math.floor(m/2);	
        double[] reF = new double[r+m-1];	
        int pin =0;
        for (int i=1-m2; i<r+m2;i++) { reF[pin] = i; pin++; }
        double[] res = Reflect(reF, 0.05, r+0.05);	
		for (int i=0;i<res.length;i++) res[i] = (int)Math.ceil(res[i]);
        double[] res2 = new double[res.length];
        for (int i=0;i<res2.length;i++) { int pin2 = (int) res[i]-1; res2[i] = data[pin2]; }	
        double[] c = convolve(res2, filter, new double[res2.length-filter.length+1]);
     return c;	
    }
	
    public static double[] coldFilt(double[] X, double[] ha, double[] hb) {		
        int r = X.length; int m = ha.length; int m2 = (int) Math.floor(m/2);		
        double[] reF = new double[r+m+m+1];
        int pin =0;
        for (int i=1-m; i<r+m+1;i++) { reF[pin] = i; pin++; }		
        double[] reF2 = new double[reF.length-1];
        for (int i=0;i<reF2.length;i++) reF2[i] = reF[i];	
        double[] res = Reflect(reF2, 0.05, r+0.05);	
        for (int i=0;i<res.length;i++) res[i] = (int)Math.ceil(res[i]);	
        double[] hao = new double[m2]; double[] hae = new double[m2]; double[] hbo = new double[m2]; double[] hbe = new double[m2];		
        int pin1 =0; int pin2=1;		
        for (int i=0;i<m2;i++) {
            hao[i] = ha[pin1]; hae[i] = ha[pin2]; hbo[i] = hb[pin1]; hbe[i] = hb[pin2];
            pin1+=2; pin2+=2;
        }		
        Vector t = new Vector();
        int pin3 = 6;
        for (int i=6;pin3<r+2*m-2;i++) { t.add(pin3); pin3 = pin3 +4; }
        t.add(pin3);
		int r2 = r/2; double[] Y = new double[r2];	
        double sum =0;	
        for (int i=0;i<ha.length;i++)  sum+=ha[i]*hb[i];		
        double[] s1 = new double[r2/2]; double[] s2 = new double[r2/2];
        pin1 = 1;	
        if (sum>0) {	
            for (int i=0;i<r2/2;i++) { s1[i] = pin1; s2[i] = pin1+1; pin1+=2; }	
        } else {
            for (int i=0;i<r2/2;i++) { s2[i] = pin1; s1[i] = pin1+1; pin1+=2; }
        }		
        double[] Fin1 = new double[t.size()]; double[] Fin2 = new double[t.size()]; double[] Fin3 = new double[t.size()]; double[] Fin4 = new double[t.size()];
        for (int i=0;i<t.size();i++) {			
            Object p= t.get(i); int pp = Integer.parseInt("" + p);            	
            int ppp = (int) res[pp-2]; if (ppp==X.length) { ppp = X.length-1; } if (ppp<0) { ppp=0; }
            Fin1[i] = X[ppp-1];
            ppp = (int) res[pp-4]; if (ppp==X.length) { ppp = X.length-1; } if (ppp<0) { ppp=0; }
            Fin2[i] = X[ppp-1];
            ppp = (int) res[pp-1]; if (ppp==X.length) { ppp = X.length-1; } if (ppp<0) { ppp=0; }
            Fin3[i] = X[ppp-1];
            ppp = (int) res[pp-3]; if (ppp==X.length) { ppp = X.length-1; } if (ppp<0) { ppp=0; }
            Fin4[i] = X[ppp-1];	
        }	
        double[] temp1 =  convolve(Fin1,hao,new double[Fin1.length-hao.length+1]); double[] temp2 =  convolve(Fin2,hae,new double[Fin2.length-hae.length+1]);
        double[] temp3 =  convolve(Fin3,hbo,new double[Fin3.length-hbo.length+1]); double[] temp4 =  convolve(Fin4,hbe,new double[Fin4.length-hbe.length+1]);				
        for (int i=0;i<r2/2;i++) {
            int p1 = (int) s2[i]-1; int p2 = (int) s1[i]-1;
            Y[p1] = temp1[i] + temp2[i]; Y[p2] = temp3[i] + temp4[i];       
        }	
    return Y;		
    }
	
    public static double[] coliFilt(double[] X, double[] ha, double[] hb) {
        int r = X.length; int m = ha.length; int m2 = (int) Math.floor(m/2); double[] Y = new double[r*2];	
        double[] hao = new double[m2]; double[] hae = new double[m2];
        double[] hbo = new double[m2]; double[] hbe = new double[m2];	
        int pin1 =0; int pin2=1;		
        for (int i=0;i<m2;i++) {
            hao[i] = ha[pin1]; hae[i] = ha[pin2];		
            hbo[i] = hb[pin1]; hbe[i] = hb[pin2];
            pin1+=2; pin2+=2;
        }				
        double[] reF = new double[r+m2+m2+1];	
        int pin =0;
        for (int i=1-m2; i<r+m2+1;i++) { reF[pin] = i; pin++;  }		
        double[] reF2 = new double[reF.length-1];	
        for (int i=0;i<reF2.length;i++) reF2[i] = reF[i];	
        double[] res = Reflect(reF2, 0.05, r+0.05);
        for (int i=0;i<res.length;i++) res[i] = (int)Math.ceil(res[i]);		
        double[] res2 = new double[res.length];	
        for (int i=0;i<res2.length;i++) { int pin3 = (int) res[i]-1; res2[i] = X[pin3]; }	
        if (m2 % 2 == 0) {		
			Vector t = new Vector(); int ppin = 4;
			for (int i=4;i<m+r/2-4;i++) { t.add(ppin); ppin=ppin+2; }		
			double sum =0;
			for (int i=0;i<ha.length;i++) sum+=ha[i]*hb[i];
			double[] ta = new double[t.size()]; double[] tb = new double[t.size()];
			
			if (sum>0) {	
				for (int i=0;i<t.size();i++) { ta[i] = Double.parseDouble("" + t.get(i)); tb[i] = Double.parseDouble("" + t.get(i))-1; }			
			} else {
                for (int i=0;i<t.size();i++) { ta[i] = Double.parseDouble("" + t.get(i))-1; tb[i] = Double.parseDouble("" + t.get(i)); }
			}		
			double[] t1 = new double[ta.length]; double[] t2 = new double[ta.length];			
			double[] t3 = new double[ta.length]; double[] t4 = new double[ta.length];			
			for (int i=0;i<ta.length;i++) {	
				int p1 = (int) ta[i]; int p2 = (int) tb[i];		
				t1[i] = res2[p1-1]; t2[i] = res2[p2-1];	
				int p3 = (int) ta[i]-2; int p4 = (int) tb[i]-2;	
				t3[i] = res2[p2-1]; t4[i] = res2[p3-1];	
			}		
			double[] temp1 =  convolve(t4,hae,new double[t1.length-hbo.length+1]);
			double[] temp2 =  convolve(t3,hbe,new double[t1.length-hbo.length+1]);
			double[] temp3 =  convolve(t2,hao,new double[t1.length-hbo.length+1]);
			double[] temp4 =  convolve(t1,hbo,new double[t1.length-hbo.length+1]);	
			int mPin = 0;
			for (int i=0;i<temp1.length;i++) {
				Y[mPin] = temp1[i]; Y[mPin+1] = temp2[i]; Y[mPin+2] = temp3[i]; Y[mPin+3] = temp4[i]; 
				mPin = mPin+4;
			}	
        } else {			
			Vector t = new Vector(); int ppin = 3;
			for (int i=4;i<m+r/2-5;i++) { t.add(ppin); ppin=ppin+2; }		
			double sum =0;
			for (int i=0;i<ha.length;i++)  sum+=ha[i]*hb[i];
			double[] ta = new double[t.size()]; double[] tb = new double[t.size()];		
			if (sum>0) {		
				for (int i=0;i<t.size();i++) { ta[i] = Double.parseDouble("" + t.get(i)); tb[i] = Double.parseDouble("" + t.get(i))-1; }		
			} else {
                for (int i=0;i<t.size();i++) { ta[i] = Double.parseDouble("" + t.get(i))-1; tb[i] = Double.parseDouble("" + t.get(i)); }
			}			
			double[] t1 = new double[ta.length]; double[] t2 = new double[ta.length];	
			double[] t3 = new double[ta.length]; double[] t4 = new double[ta.length];	
			for (int i=0;i<ta.length;i++) {				
				int p1 = (int) ta[i]; int p2 = (int) tb[i];			
				t1[i] = res2[p1-1]; t2[i] = res2[p2-1];	
				int p3 = (int) ta[i]-2; int p4 = (int) tb[i]-2;		
				t3[i] = res2[p2-1]; t4[i] = res2[p3-1];
			}		
			double[] temp1 =  convolve(t2,hao,new double[t1.length-hbo.length+1]);
			double[] temp2 =  convolve(t1,hbo,new double[t1.length-hbo.length+1]);
			double[] temp3 =  convolve(t2,hae,new double[t1.length-hbo.length+1]);
			double[] temp4 =  convolve(t1,hbe,new double[t1.length-hbo.length+1]);			
			int mPin = 0;
			for (int i=0;i<temp1.length;i++) {
				Y[mPin] = temp1[i]; Y[mPin+1] = temp2[i]; Y[mPin+2] = temp3[i]; Y[mPin+3] = temp4[i];
				mPin = mPin+4;
			}	
        }	
    return Y;	
    }
	
    
   public static void main(String[] args) {	
		String file ="";		
		for (int u=0;u<args.length; u++) {	
			int pin = u+1; String s = args[u]; 		
			if (s.equals("-f")) file = args[pin];
			if (s.equals("-fa")) factor = Double.parseDouble(""+args[pin]);
			if (s.equals("-si")) sized = Double.parseDouble(""+args[pin]);
		}	
		ArrayList[] data = Input(file); double[] ddd = new double[data[0].size()];   	
    	for(int i=0;i<data[0].size();i++) ddd[i] = Double.parseDouble(""+data[3].get(i));	
		double[] wav = dtWaveXfm(ddd, 1); double[] ws = rtab(wav,ddd);
		for (int i=0;i<ddd.length;i++) System.out.println(ddd[i] + "\t" + wav[i] + "\t" + ws[i]);
	}
}






