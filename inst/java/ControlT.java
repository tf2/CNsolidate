
import java.util.ArrayList;
import java.util.Vector;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Vector;

/**
 *
 * @author tf2
 */
public class ControlT {

	public String file="";
	public static double minRatioValue = 0.35;
	public static double tr1 = 2; public static int len1 = 25; public static int len2 = 5; public static int wlen1 = 100;
    public static int wlen2 = 500; public static double tr2 = 2; public static int wlen3 = 10; public static int siZ = 500;
    
    public ControlT(String filename) {
    	this.file = filename; int index = file.lastIndexOf(".");
		String inName = file;String ouName = file.substring(0, index) + "_walk.temp";
		Handles hand = new Handles(file, ouName); ArrayList[] data = hand.getData(4);
		int[] chrs = new int[data[0].size()]; int[] starts = new int[data[0].size()];
		int[] stops = new int[data[0].size()]; double[] ddd = new double[data[0].size()];

		for(int i=0;i<chrs.length;i++) {
			chrs[i] = Integer.parseInt( "" + data[0].get(i)); starts[i] = Integer.parseInt( "" + data[1].get(i));
			stops[i] = Integer.parseInt( "" + data[2].get(i)); ddd[i] = Double.parseDouble( "" + data[3].get(i));
		}
        double[] r1 = RunItApp(ddd); r1 = exMath.runMedian(r1, 1001);                 
        double[] r3 = exMath.runMedianR(r1, 1001); r3 =getRunningLower(r3, 0.68);
        double[] r4 = findMedians(ddd, r3); r4 = checkSmallCall(ddd,r4, r3);

        r4 = findMedians(ddd, r4);   
        for(int i=0;i<r4.length;i++) { if (Math.abs(r4[i]) < minRatioValue) r4[i] = 0; } 
        hand.Out(chrs, starts, stops, ddd, r3, r4, ouName);
    }

    public static double[] findInterval(double[] l) {    
        double[] v = new double[l.length]; double d = exMath.dLRs(l);
        double t = d*tr1; int e = v.length-len1;
        for (int x=0;x<e;x++) {
        	int s =0; v[x] =0; double m = Math.abs(l[x]);
            if(m>t) {
            	for(int i=0;i<len1;i++) {
                	double m2 = Math.abs(l[x+i]);
                    if (m2>t) s=s+1;
                }
            } if (s > len2) v[x] = 1;
        }
        return v;
    }

     public static double[] mergeInterval(double[] v) {
         double[] v1 = new double[v.length];double[] v2 = new double[v.length];
         double[] vv = new double[v.length]; int e = v.length-wlen1;
         for (int x=0;x<e;x++) {
             int s =0;
             for (int i=0;i<wlen1;i++) { if (v[i+x] == 1) s=1; }
             v1[x] =s;
         }
         for (int x=e;x>=wlen1;x--) {
             int s =0;
             for (int i=wlen1-1;i>=0;i--) { if (v[i+x] == 1) s=1; }
             v2[x] =s;
         }
         for (int i=0;i<vv.length;i++) {
             vv[i] = v[i]; int t = (int) (v1[i] + v2[i]);
             if(t == 1) vv[i]=0;
             else if (t == 2) vv[i]=1;
         }
         for (int k=vv.length-wlen1;k>=0;k--) {
            if (vv[k]==0) { for (int j=k;j<k+wlen1;j++) vv[j]=0; }
        }
         return vv;
     }

     public static double[] filterInterval(double[] v) {
         double[] cc = new double[v.length]; int e = v.length-wlen2;
         for (int x=wlen2;x<e;x++) {
             int s1 =1; int s2 =1; cc[x] = v[x];   
                 for(int i=x-wlen2;i<x;i++) { if (v[i] == 0) s1 = 0; }
                 for (int k=x+wlen2;k>=x;k--) { if (v[k] == 0) s2 = 0; }
             if(s1 == 0 || s2 == 0)  cc[x] =0;
         }
        for (int k=0;k<cc.length;k++) { if (cc[k]==1) { for (int j=k;j>=k-wlen2;j--) cc[j] =1; } }
        for (int k=cc.length-1;k>=0;k--) { if (cc[k]==1) { for (int j=k;j<k+wlen2;j++) cc[j] =1; } }
        return cc;
     }

     public static double[] findMedians(double[] l, double[] v) {
         boolean t = true; int k=0;int pin1 =0;int pin2 = 0;
         Vector meds = new Vector();double[] r = new double[l.length]; v[0] = 0;
         do {
         	while(pin1 < v.length-1 && v[pin1] == 0) pin1++;
            	if (v[pin1]==1) {
                    double[] me = new double[pin1-pin2]; System.arraycopy(l,k , me, 0, me.length);
                    double md = exMath.median(me); meds.add(md);
                    k = pin1; pin2 = pin1;
                }
            while(pin1 < v.length-1 && v[pin1] == 1) pin1++;           
                if (v[pin1]==0) {
                    double[] me = new double[pin1-pin2]; System.arraycopy(l,k, me, 0, me.length);
                    double md = exMath.median(me); meds.add(md);
                    k=pin1;pin2 = pin1;
                }
           if (pin1 >= v.length-1) {  t = false; break; }        
         } while(t);

         t = true;
         k=0; pin1 =0; pin2 = 0;
         ArrayList[] MEDS = new ArrayList[meds.size()];
         for (int i=0;i<MEDS.length;i++) MEDS[i] = new ArrayList();
            do {
                if (pin1 >= v.length-1) {  t = false; break; }
                while(pin1 < v.length-1 && v[pin1] == 0) pin1++;
                	if (v[pin1]==1) {
                    	for (int j=k;j<pin1;j++) {
                        r[j] = Double.parseDouble("" + meds.get(pin2)); MEDS[pin2].add(l[j]);
                    } k=pin1;pin2++;
                }
                while(pin1 < v.length-1 && v[pin1] == 1) pin1++;
                if (v[pin1]==0) {
                    for (int j=k;j<pin1;j++) {
                        r[j] = Double.parseDouble("" + meds.get(pin2)); MEDS[pin2].add(l[j]);
                    } k=pin1; pin2++;
                }           
            } while(t);

            r[r.length-1] = Double.parseDouble("" + meds.get(pin2-1)); MEDS[pin2-1].add(l[l.length-1]);
            r = extendFilter(l, r, 20);
         return r;
     }

     public static double[] findMeans(double[] l, double[] v) {
         boolean t = true; int k=0; int pin1 =0;int pin2 = 0;
         Vector meds = new Vector(); double[] r = new double[l.length]; v[0] = 0;
         do {
         	while(pin1 < v.length-1 && v[pin1] == 0) pin1++;
                if (v[pin1]==1) {
                    double[] me = new double[pin1-pin2]; System.arraycopy(l,k , me, 0, me.length);
                    double md = exMath.mean(me); meds.add(md);
                    k = pin1; pin2 = pin1;
                }
            while(pin1 < v.length-1 && v[pin1] == 1) pin1++;
                if (v[pin1]==0) {
                    double[] me = new double[pin1-pin2]; System.arraycopy(l,k, me, 0, me.length);
                    double md = exMath.mean(me); meds.add(md);
                    k=pin1; pin2 = pin1;
                }
            if (pin1 >= v.length-1) {  t = false; break; }
         } while(t);

         t = true;
         k=0; pin1 =0; pin2 = 0;
         ArrayList[] MEDS = new ArrayList[meds.size()];
         for (int i=0;i<MEDS.length;i++) MEDS[i] = new ArrayList();
           	do {
                if (pin1 >= v.length-1) {  t = false; break; }
                while(pin1 < v.length-1 && v[pin1] == 0) pin1++;
                if (v[pin1]==1) {
                    for (int j=k;j<pin1;j++) {
                        r[j] = Double.parseDouble("" + meds.get(pin2)); MEDS[pin2].add(l[j]);
                    } k=pin1; pin2++;
                }
                while(pin1 < v.length-1 && v[pin1] == 1) pin1++;
                if (v[pin1]==0) {
                    for (int j=k;j<pin1;j++) {
                        r[j] = Double.parseDouble("" + meds.get(pin2)); MEDS[pin2].add(l[j]);
                    } k=pin1; pin2++;
                }
            } while(t);

            r[r.length-1] = Double.parseDouble("" + meds.get(pin2-1)); MEDS[pin2-1].add(l[l.length-1]);
            r = extendFilter(l, r, 20);
         return r;
     }


     public static double[] extendFilter(double[] l, double[] r, int wen) {
     	for (int i=0;i<l.length-wen;i++) {
        	boolean ch = true;
            if (r[i]<0) { for(int j=i;j<i+wen;j++) { if(l[j]<r[i]) ch = false; }
            } else if (r[i]>0) { for(int j=i;j<i+wen;j++) { if(l[j]>r[i]) ch = false; }
            }
            if (ch) {
            	double p = 2;
                for(int kk=i;kk<r.length;kk++) { if (r[i]!=r[kk]) { p = r[kk]; break; } }
                	r[i] = p;
                }
            }
        for (int i=l.length-1;i>wen;i--) {
        	boolean ch = true;
           	if (r[i]<0) { for(int j=i;j>i-wen;j--) { if(l[j]<r[i]) ch = false; }
            } else if (r[i]>0) { for(int j=i;j>i-wen;j--) { if(l[j]>r[i]) ch = false; }
            }
            if (ch) {
            	double p = 2;
                for(int kk=i;kk>0;kk--) { if (r[i]!=r[kk]) { p = r[kk]; break; } }
                r[i] = p;
                }
            }
     	return r;
     }

     public static double[] fineFilter(double[] l, double[] v) {
         boolean t = true; int k=0; int pin1 =0; int pin2 = 0;
         Vector meds = new Vector(); v[0] = 0;
         do {
         	while(pin1 < v.length-1 && v[pin1] == 0) pin1++;
                if (v[pin1]==1) {
                    double[] me = new double[pin1-pin2]; System.arraycopy(l,k , me, 0, me.length);
                    double md = exMath.median(me); meds.add(md);
                    k = pin1; pin2 = pin1;
                }
            while(pin1 < v.length-1 && v[pin1] == 1) pin1++;           
                if (v[pin1]==0) {
                    double[] me = new double[pin1-pin2]; System.arraycopy(l,k, me, 0, me.length);
                    double md = exMath.median(me); meds.add(md);
                    k=pin1;  pin2 = pin1;
                }
            if (pin1 >= v.length-1) {  t = false; break; }              
         } while(t);
         t = true;
         k=0; pin1 =0;  pin2 = 0;
         ArrayList[] MEDS = new ArrayList[meds.size()];
         for (int i=0;i<MEDS.length;i++) MEDS[i] = new ArrayList();
            do {
                if (pin1 >= v.length-1) {  t = false; break; }
                while(pin1 < v.length-1 && v[pin1] == 0) pin1++;
                if (v[pin1]==1) { 
                	for (int j=k;j<pin1;j++) { MEDS[pin2].add(l[j]); }
                    k=pin1; pin2++;
                }
                while(pin1 < v.length-1 && v[pin1] == 1) pin1++;
                if (v[pin1]==0) {
                    for (int j=k;j<pin1;j++) { MEDS[pin2].add(l[j]); }
                    k=pin1; pin2++;
                }           
            } while(t);
         MEDS[pin2-1].add(l[l.length-1]); ArrayList[] NEW_MEDS = new ArrayList[MEDS.length];
         for (int i=0;i<MEDS.length;i++) NEW_MEDS[i] = new ArrayList();
         Vector pinner = new Vector();
         for (int i=0;i<MEDS.length;i++) {
              double[] temp = new double[MEDS[i].size()];
                    for (int j=0;j<MEDS[i].size();j++) temp[j] = Double.parseDouble("" + MEDS[i].get(j));;
              double d = exMath.dLRs(temp)*tr2; double rp68 = exMath.quantile(exMath.Absol(temp), 0.68);
                    for (int a=0;a<temp.length;a++) {
                        NEW_MEDS[i].add(temp[a]); boolean check = true;
                        if (a<temp.length-wlen3) {
                            for (int tt=a+1;tt<a+wlen3;tt++) {
                                double dif = temp[a] - temp[tt];
                                if (dif<d) check = false;
                            }
                        } else if (a>temp.length-wlen3) {
                            for (int tt=temp.length-1;tt>temp.length-wlen3;tt--) {
                                double dif = temp[a] - temp[tt];
                                if (dif<d) check = false;
                            }
                        }              
                        if(check) pinner.add(NEW_MEDS[i].get(a)); else pinner.add(0);
                    }
                }
         double[] medis = new double[pinner.size()];
         for (int i=0;i<pinner.size();i++) medis[i] = Double.parseDouble("" + pinner.get(i));
         medis[0] = 0; t = true;
         k=0; pin1 =0; pin2 = 0;
         double[] r = new double[l.length]; Vector fin = new Vector();
         do {
         	while(pin1 < medis.length-1 && medis[pin1] == 0) pin1++;
                if (medis[pin1]!=0) {
                    double[] me = new double[pin1-pin2]; System.arraycopy(l,k , me, 0, me.length);
                    double md = exMath.median(me); fin.add(md);
                    k = pin1;pin2 = pin1;
                }
            while(pin1 < medis.length-1 && medis[pin1] != 0) pin1++;
                if (medis[pin1]==0) {
                    double[] me = new double[pin1-pin2]; System.arraycopy(l,k, me, 0, me.length);
                    double md = exMath.median(me); fin.add(md);
                    k=pin1; pin2 = pin1;
                }
            if (pin1 >= medis.length-1) {  t = false; break; }
         } while(t);
         t = true;
         k=0; pin1 =0;  pin2 = 0;
            do {
                if (pin1 >= medis.length-1) {  t = false; break; }
                while(pin1 < medis.length-1 && medis[pin1] == 0) pin1++;
                if (medis[pin1]!=0) {
                    for (int j=k;j<pin1;j++) r[j] = Double.parseDouble("" + fin.get(pin2));
                    k=pin1; pin2++;
                }
                while(pin1 < medis.length-1 && medis[pin1] !=0) pin1++;
                	if (medis[pin1]==0) { for (int j=k;j<pin1;j++) r[j] = Double.parseDouble("" + fin.get(pin2));
               		k=pin1; pin2++;
                }
            } while(t);
            r[r.length-1] = Double.parseDouble("" + fin.get(pin2-1));
            double[] a = exMath.Absol(r); double d = exMath.quantile(a, 0.68);
            a = exMath.Absol(l); double dl = exMath.quantile(a, 0.68);
            for (int i=0;i<r.length-1;i++) { double diff = Math.abs(r[i] - r[i+1]); if (diff > d) r[i] =0; }
            for(int i=1;i<r.length;i++) { if (r[i]==0) r[i] = r[i-1]; }
            
     	return r;
     }

  	public static double[] getRunningLower(double[] r, double per) {
         double[] d = new double[r.length];
         for(int i=0;i<r.length-1;i++) d[i] = Math.abs(r[i] - r[i+1]);
         d[d.length-1] = 0; Vector v = new Vector();
         for (int i=0;i<d.length;i++) { if(d[i] > 0) v.add(d[i]); }
         double[] t = new double[v.size()];
         for (int i=0;i<t.length;i++) t[i] = Double.parseDouble("" + v.get(i));

         double q = exMath.quantile(t, per);
         for (int i=0;i<d.length; i++) {
             if(d[i]<q) d[i] = 0; else d[i] = 1;
         }
     	return d;
     }

     public static double[] checkSmallCall(double[] l, double[] r, double[] rr) {
         double[] a = exMath.Absol(l); double t = exMath.quantile(a, 0.68)*tr1;
         for (int i=0;i<l.length;i++) { double d = Math.abs(l[i] - r[i]); if (d>t) rr[i] = 1; }
         double[] p = new double[rr.length];
         for(int i=1;i<rr.length-1;i++) {
            p[i] = rr[i];
            if(rr[i]==1) {
                if(rr[i]==rr[i-1]) p[i] = 0;
                if (rr[i+1]==0) p[i] = 1;
            } 
         }
     	return p;
     }

    public static double[][] RunMeanWalk(double[] data) {
        double[][] retData = new double[data.length][3]; int[] windows = { 5,10,20,50,100 };
        for (int i=0;i<windows.length;i++) {
            int pin =0;
            for (int j=0;j<data.length/windows[i];j++) {
                double[] inter = new double[windows[i]];System.arraycopy(data, pin, inter, 0, inter.length);
                double m = exMath.mean(inter);double v = exMath.Var(inter);   
                for (int k =0;k<windows[i];k++) {
                    retData[pin+k][0] = pin; retData[pin+k][1] = m;retData[pin+k][2] = v;
               } pin += windows[i];
            }
        }
    	return retData;
    }

    public static Object[] RandomWalk(double[] ordata, double[] data, int win_len) {     
        double[] data_filter = new double[data.length];  
        double m = exMath.mean(data);int win_over =1;
        for (int i=0;i<data.length;i++) data_filter[i] = data[i] -m;     
        double[] data_filter_int = new double[data.length]; data_filter_int[0] = data_filter[0];
        for (int i=1;i<data.length;i++) data_filter_int[i] = data_filter[i-1] + data_filter[i];   
        double S = Math.floor((data.length - win_len) /win_over+1); double[] keep = new double[data.length];        
        for (int i=1;i<S;i++) {      
            int fro = (i-1)*win_over+1; double p = (i+(win_len-1)/2); int pp = (int) p;

            double[] temp = new double[win_len]; System.arraycopy(data, fro, temp, 0, temp.length);      
            keep[pp] =exMath.Var(temp);
        }
        double[] newkeep = new double[keep.length]; double[] result = new double[keep.length]; double[] newresult = new double[keep.length];
        for (int i=0;i<keep.length;i++) newkeep[i] = keep[i];
        int st = (win_len+1)/2; int sto = keep.length-(win_len-1)/2;
        for (int i=st;i<sto;i++) {
            double[] v1 = new double[1]; double[] v2 = new double[1]; int pin = i-1;
            System.arraycopy(data, pin, v1, 0, v1.length); System.arraycopy(data, i, v2, 0, v2.length);
            if (exMath.Var(v1) >= exMath.Var(v2)) newkeep[i] = data_filter[i+1];
            else newkeep[i] = data_filter[i-1];
            result[i] = newkeep[i]*keep[i]; newresult[i] = newkeep[i]*keep[i];
        }
        double thres_set = 0.05;
        for (int i=0;i<result.length;i++) { if (Math.abs(result[i])>thres_set) { newresult[i] = newkeep[i]; } }
        double cut_mean = exMath.mean(newresult); Object[] fin = new Object[newresult.length];
        fin[0] = "NaN"; fin[1] = "NaN"; fin[2] = "NaN";
        for (int i=2;i<newresult.length;i++) {
            int pin1 = i-1; int pin2 = i-2;
            if (i>0 && Math.abs(newresult[i]) > Math.abs(cut_mean)*5 && Math.abs(newresult[pin1]) >  Math.abs(cut_mean)*5 && Math.abs(newresult[pin2]) >  Math.abs(cut_mean)*5) {
               fin[i] = newresult[i]; fin[pin1] = newresult[pin1]; fin[pin2] = newresult[pin2];
            } else { fin[i] = "NaN"; }
        }
        for (int i=1;i<fin.length;i++) fin[i-1] = fin[i];
    	return fin;
    }

    public static double[] RunItApp(double[] l) {
         double[] rr = new double[l.length];
         double[] v = findInterval(l); double[] vv = mergeInterval(v);
         double[] vvv = filterInterval(vv); double[] r = findMedians(l, vvv);
         rr = fineFilter(l, vvv);
         Object[] v1 = new Object[r.length]; Object[] v2 = new Object[r.length];
         for ( int i=0;i<r.length;i++) { v1[i] = r[i]; v2[i] = rr[i]; }
        return rr;
    }

public static void main(String[] args) {
		String file=""; String iDir = ""; String oDir = "";	
		for (int u=0;u<args.length; u++) {			
			int pin = u+1; String s = args[u]; 		
			if (s.equals("-f")) file = args[pin];			
			if (s.equals("-mr")) minRatioValue = Double.parseDouble("" + args[pin]);		
			if (s.equals("-i")) iDir = args[pin];	
			if (s.equals("-o")) oDir = args[pin];		
			if( s.equals("-t1") ) tr1 = Integer.parseInt("" + args[pin]);
			if( s.equals("-t2") ) tr2 = Integer.parseInt("" + args[pin]);
			if( s.equals("-l1") ) len1 = Integer.parseInt("" + args[pin]);		
			if( s.equals("-l2") ) len2 = Integer.parseInt("" + args[pin]);		
			if( s.equals("-w1") ) wlen1 = Integer.parseInt("" + args[pin]);			
			if( s.equals("-w2") ) wlen2 = Integer.parseInt("" + args[pin]);    		
			if( s.equals("-w3") ) wlen3 = Integer.parseInt("" + args[pin]);	
			if( s.equals("-s") ) siZ = Integer.parseInt("" + args[pin]);			
		}
		ControlT conT = new ControlT(file);
    }

}
