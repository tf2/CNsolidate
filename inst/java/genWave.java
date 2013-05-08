import java.io.*;
import java.util.*;
import java.lang.Math.*;
import java.util.ArrayList;


public class genWave {

	/******************************************************************************/
	/******************************** Wavelets ************************************/
	
	final static double[] g0a = {-0.0046, -0.0054, 0.0170, 0.0238, -0.1067, 0.0119, 0.5688,
	0.7561, 0.2753, -0.1172, -0.0389, 0.0347, -0.0039, 0.0033};
	
    final static double[] g0b = {0.0033, -0.0039, 0.0347, -0.0389, -0.1172, 0.2753, 0.7561,
	0.5688, 0.0119, -0.1067, 0.0238, 0.0170, -0.0054, -0.0046};
	
    final static double[] g1a = {-0.0033, -0.0039, -0.0347, -0.0389, 0.1172, 0.2753, -0.7561,
	0.5688, -0.0119, -0.1067, -0.0238, 0.0170, 0.0054, -0.0046};
	
    final static double[] g1b = {-0.0046, 0.0054, 0.0170, -0.0238, -0.1067, -0.0119, 0.5688,
	-0.7561, 0.2753, 0.1172, -0.0389, -0.0347, -0.0039, -0.0033};
	
    final static double[] g0o = {0.0001, 0, -0.0013, -0.0019, 0.0072, 0.0239, -0.0556, -0.0517, 0.2998, 0.5594,
	0.2998, -0.0517, -0.0556, 0.0239, 0.0072, -0.0019, -0.0013, 0, 0.0001};
	
    final static double[] g1o = {-0.0018, 0, 0.0223, 0.0469, -0.0482, -0.2969,
	0.5555, -0.2969, -0.0482, 0.0469, 0.0223, 0, -0.0018};
	
	
	
    final static double[] h0a = {0.0033, -0.0039, 0.0347, -0.0389, -0.1172, 0.2753, 0.7561,
	0.5688, 0.0119, -0.1067, 0.0238, 0.0170, -0.0054, -0.0046};
	
    final static double[] h0b = {-0.0046, -0.0054, 0.0170, 0.0238, -0.1067, 0.0119, 0.5688,
	0.7561, 0.2753, -0.1172, -0.0389, 0.0347, -0.0039, 0.0033};
	
    final static double[] h1a = {-0.0046, 0.0054, 0.0170, -0.0238, -0.1067, -0.0119, 0.5688,
	-0.7561, 0.2753, 0.1172, -0.0389, -0.0347, -0.0039, -0.0033};
	
    final static double[] h1b = {-0.0033, -0.0039, -0.0347, -0.0389, 0.1172, 0.2753, -0.7561,
	0.5688, -0.0119, -0.1067, -0.0238, 0.0170, 0.0054, -0.0046};
	
	final static double[] h0o = {-0.0018, 0, 0.0223, -0.0469, -0.0482, 0.2969,
	0.5555, 0.2969, -0.0482, -0.0469, 0.0223, 0, -0.0018};
	
    final static double[] h1o = {-0.0001, 0, 0.0013, -0.0019, -0.0072, 0.0239, 0.0556, -0.0517, -0.2998, 0.5594,
	-0.2998, -0.0517, 0.0556, 0.0239, -0.0072, -0.0019, 0.0013, 0, -0.0001};
	
	
	
	public genWave() {
	
	}
	
	
	public static double mean(Object[] nn) {
		
		double[] m = new double[nn.length];
		for(int i=0;i<m.length;i++) {
			m[i]=Double.parseDouble(""+nn[i]);
		}
		
		double men = 0;
		int l = m.length;
		
		for (int i = 0 ; i<m.length ; i++) {
			men = men + m[i];     
		}
		
		return men / l;
	}
	
	public static double[] Absol(double[] values) {
        
        double[] a = new double[values.length];
		
        for (int i =0; i<values.length ; i++) {
            a[i] = Math.abs(values[i]);
        }
        
        return a;    
    }
	
	public static double quantile(double[] values, double quantile) {
        if(values == null)
            throw new IllegalArgumentException("Values cannot be null.");
        if(quantile < 0.0 || quantile > 1.0)
            throw new IllegalArgumentException("Quantile must be between 0.0 and 1.0");
		
        double[] copy = new double[values.length];
        System.arraycopy(values, 0, copy, 0, copy.length);
        java.util.Arrays.sort(copy);
        int index = (int) (copy.length * quantile);
        return copy[index];
    }
	
	/******************************************************************************/
	/********************************** gWave *************************************/
	
	
    public static double[] dtWaveXfm(double[] data, int nLevels) {
		
        int L = data.length;
		boolean odd = false;
		
        if (L % 2 == 1 ) {
			odd=true;
			L++;
        }
		
		double[] data2 = new double[L];
		for(int i=0;i<L-1;i++) {
			data2[i] = data[i];
		}
		if(odd) {
			data2[L-1] = 0;
		}
		
        double[] Hi = colFilter(data2, h1o);
        double[] Lo = colFilter(data2, h0o);
        
        for (int i=0;i<nLevels;i++) {
			
            if (Lo.length<10) { break; }
			
            if (Lo.length % 4 != 0 ) {
                Lo = extend(Lo);
            }           
			
            Hi = coldFilt(Lo,h1b,h1a);
            Lo = coldFilt(Lo,h0b,h0a);
			
			Lo = dtWaveIfm(Lo, Hi, nLevels);
			
        }          
		
        return Lo;
		
    }
	
	
    public static double[] dtWaveIfm(double[] Yl, double[] Yh, int nLevels) {
		
		double[] Lo = Yl;
		
		for (int i =0;i<nLevels;i++) {
			Lo = coliFilt(Yl, g0b, g0a);
		}
		
        return Lo;
		
    }
	
	
    public static double[] colFilter(double[] data, double[] filter) {
		
        int r = data.length;
        int m = filter.length;
        int m2 = (int) Math.floor(m/2);
		
        double[] reF = new double[r+m-1];
		
        int pin =0;
        for (int i=1-m2; i<r+m2;i++) {
            reF[pin] = i;
            pin++;
        }
		
        double[] res = Reflect(reF, 0.05, r+0.05);
		
		for (int i=0;i<res.length;i++) {
			res[i] = (int)Math.ceil(res[i]);
		}
		
        double[] res2 = new double[res.length];
		
        for (int i=0;i<res2.length;i++) {
            int pin2 = (int) res[i]-1;
            res2[i] = data[pin2];
        }
		
        double[] c = convolve(res2, filter, new double[res2.length-filter.length+1]);
		
		
        return c;
		
    }
	
    public static double[] coldFilt(double[] X, double[] ha, double[] hb) {
		
        int r = X.length;
        int m = ha.length;
        int m2 = (int) Math.floor(m/2);
		
        double[] reF = new double[r+m+m+1];
		
        int pin =0;
        for (int i=1-m; i<r+m+1;i++) {
            reF[pin] = i;
            pin++;
        }
		
        double[] reF2 = new double[reF.length-1];
		
        for (int i=0;i<reF2.length;i++) {
            reF2[i] = reF[i];
        }
		
        double[] res = Reflect(reF2, 0.05, r+0.05);
		
        for (int i=0;i<res.length;i++) {
			res[i] = (int)Math.ceil(res[i]);
		}
		
        double[] hao = new double[m2];
        double[] hae = new double[m2];
        double[] hbo = new double[m2];
        double[] hbe = new double[m2];
		
        int pin1 =0;
        int pin2=1;
		
        for (int i=0;i<m2;i++) {
            hao[i] = ha[pin1];
            hae[i] = ha[pin2];
			
            hbo[i] = hb[pin1];
            hbe[i] = hb[pin2];
            pin1+=2;
            pin2+=2;
        }
		
        Vector t = new Vector();
        int pin3 = 6;
        for (int i=6;pin3<r+2*m-2;i++) {
            t.add(pin3);
            pin3 = pin3 +4;
        }
        t.add(pin3);
		
		int r2 = r/2;
		double[] Y = new double[r2];
		
		
        double sum =0;
		
        for (int i=0;i<ha.length;i++) {
            sum+=ha[i]*hb[i];
        }
		
        double[] s1 = new double[r2/2];
        double[] s2 = new double[r2/2];
        pin1 = 1;
		
        if (sum>0) {
			
            for (int i=0;i<r2/2;i++) {
                s1[i] = pin1;
                s2[i] = pin1+1;
                pin1+=2;
            }
			
        } else {
            for (int i=0;i<r2/2;i++) {
                s2[i] = pin1;
                s1[i] = pin1+1;
                pin1+=2;
            }
        }
		
        double[] Fin1 = new double[t.size()];
        double[] Fin2 = new double[t.size()];
        double[] Fin3 = new double[t.size()];
        double[] Fin4 = new double[t.size()];
		
        for (int i=0;i<t.size();i++) {
			
            Object p= t.get(i);
            int pp = Integer.parseInt("" + p);            
			
            int ppp = (int) res[pp-2]; if (ppp==X.length) { ppp = X.length-1; } if (ppp<0) { ppp=0; }
            Fin1[i] = X[ppp-1];
            ppp = (int) res[pp-4]; if (ppp==X.length) { ppp = X.length-1; } if (ppp<0) { ppp=0; }
            Fin2[i] = X[ppp-1];
            ppp = (int) res[pp-1]; if (ppp==X.length) { ppp = X.length-1; } if (ppp<0) { ppp=0; }
            Fin3[i] = X[ppp-1];
            ppp = (int) res[pp-3]; if (ppp==X.length) { ppp = X.length-1; } if (ppp<0) { ppp=0; }
            Fin4[i] = X[ppp-1];
			
        }
		
        double[] temp1 =  convolve(Fin1,hao,new double[Fin1.length-hao.length+1]);
        double[] temp2 =  convolve(Fin2,hae,new double[Fin2.length-hae.length+1]);
        double[] temp3 =  convolve(Fin3,hbo,new double[Fin3.length-hbo.length+1]);
        double[] temp4 =  convolve(Fin4,hbe,new double[Fin4.length-hbe.length+1]);
		
		
        for (int i=0;i<r2/2;i++) {
            int p1 = (int) s2[i]-1;
            int p2 = (int) s1[i]-1;
            Y[p1] = temp1[i] + temp2[i];
            Y[p2] = temp3[i] + temp4[i];       
        }
		
        return Y;
		
    }
	
    public static double[] coliFilt(double[] X, double[] ha, double[] hb) {
		
		
        int r = X.length;
        int m = ha.length;
        int m2 = (int) Math.floor(m/2);
		
        double[] Y = new double[r*2];
		
        double[] hao = new double[m2];
        double[] hae = new double[m2];
        double[] hbo = new double[m2];
        double[] hbe = new double[m2];
		
        int pin1 =0;
        int pin2=1;
		
        for (int i=0;i<m2;i++) {
            hao[i] = ha[pin1];
            hae[i] = ha[pin2];
			
            hbo[i] = hb[pin1];
            hbe[i] = hb[pin2];
            pin1+=2;
            pin2+=2;
        }
		
		
        double[] reF = new double[r+m2+m2+1];
		
        int pin =0;
        for (int i=1-m2; i<r+m2+1;i++) {
            reF[pin] = i;
            pin++;
        }
		
        double[] reF2 = new double[reF.length-1];
		
        for (int i=0;i<reF2.length;i++) {
            reF2[i] = reF[i];
        }
		
        double[] res = Reflect(reF2, 0.05, r+0.05);
		
        for (int i=0;i<res.length;i++) {
			res[i] = (int)Math.ceil(res[i]);
		}
		
        double[] res2 = new double[res.length];
		
        for (int i=0;i<res2.length;i++) {
            int pin3 = (int) res[i]-1;
            res2[i] = X[pin3];
        }
		
        if (m2 % 2 == 0) {
			
			Vector t = new Vector();
			int ppin = 4;
			for (int i=4;i<m+r/2-4;i++) {
				t.add(ppin);
				ppin=ppin+2;
			}
			
			double sum =0;
			for (int i=0;i<ha.length;i++) {
				sum+=ha[i]*hb[i];
			}
			double[] ta = new double[t.size()];
			double[] tb = new double[t.size()];
			
			if (sum>0) {
				
				for (int i=0;i<t.size();i++) {
					ta[i] = Double.parseDouble("" + t.get(i));
					tb[i] = Double.parseDouble("" + t.get(i))-1;
				}
				
			} else {
                for (int i=0;i<t.size();i++) {
					ta[i] = Double.parseDouble("" + t.get(i))-1;
					tb[i] = Double.parseDouble("" + t.get(i));
				}
			}
			
			double[] t1 = new double[ta.length];
			double[] t2 = new double[ta.length];
			
			double[] t3 = new double[ta.length];
			double[] t4 = new double[ta.length];
			
			for (int i=0;i<ta.length;i++) {
				
				int p1 = (int) ta[i];
				int p2 = (int) tb[i];
				
				t1[i] = res2[p1-1];
				t2[i] = res2[p2-1];
				
				int p3 = (int) ta[i]-2;
				int p4 = (int) tb[i]-2;
				
				t3[i] = res2[p2-1];
				t4[i] = res2[p3-1];
				
			}
			
			double[] temp1 =  convolve(t4,hae,new double[t1.length-hbo.length+1]);
			double[] temp2 =  convolve(t3,hbe,new double[t1.length-hbo.length+1]);
			double[] temp3 =  convolve(t2,hao,new double[t1.length-hbo.length+1]);
			double[] temp4 =  convolve(t1,hbo,new double[t1.length-hbo.length+1]);
			
			int mPin = 0;
			for (int i=0;i<temp1.length;i++) {
				Y[mPin] = temp1[i];
				Y[mPin+1] = temp2[i];
				Y[mPin+2] = temp3[i];
				Y[mPin+3] = temp4[i];
				mPin = mPin+4;
			}
			
        } else {
			
			Vector t = new Vector();
			int ppin = 3;
			for (int i=4;i<m+r/2-5;i++) {
				t.add(ppin);
				ppin=ppin+2;
			}
			
			double sum =0;
			for (int i=0;i<ha.length;i++) {
				sum+=ha[i]*hb[i];
			}
			double[] ta = new double[t.size()];
			double[] tb = new double[t.size()];
			
			if (sum>0) {
				
				for (int i=0;i<t.size();i++) {
					ta[i] = Double.parseDouble("" + t.get(i));
					tb[i] = Double.parseDouble("" + t.get(i))-1;
				}
				
			} else {
                for (int i=0;i<t.size();i++) {
					ta[i] = Double.parseDouble("" + t.get(i))-1;
					tb[i] = Double.parseDouble("" + t.get(i));
				}
			}
			
			double[] t1 = new double[ta.length];
			double[] t2 = new double[ta.length];
			
			double[] t3 = new double[ta.length];
			double[] t4 = new double[ta.length];
			
			for (int i=0;i<ta.length;i++) {
				
				int p1 = (int) ta[i];
				int p2 = (int) tb[i];
				
				t1[i] = res2[p1-1];
				t2[i] = res2[p2-1];
				
				int p3 = (int) ta[i]-2;
				int p4 = (int) tb[i]-2;
				
				t3[i] = res2[p2-1];
				t4[i] = res2[p3-1];
			}
			
			double[] temp1 =  convolve(t2,hao,new double[t1.length-hbo.length+1]);
			double[] temp2 =  convolve(t1,hbo,new double[t1.length-hbo.length+1]);
			double[] temp3 =  convolve(t2,hae,new double[t1.length-hbo.length+1]);
			double[] temp4 =  convolve(t1,hbe,new double[t1.length-hbo.length+1]);
			
			int mPin = 0;
			for (int i=0;i<temp1.length;i++) {
				Y[mPin] = temp1[i];
				Y[mPin+1] = temp2[i];
				Y[mPin+2] = temp3[i];
				Y[mPin+3] = temp4[i];
				mPin = mPin+4;
			}
			
        }
		
        return Y;
		
    }
	
	
    public static double[] Reflect(double[] dataX, double minx, double maxx) {
		
        double[] dataY = new double[dataX.length];
        
        for (int i=0;i<dataX.length;i++) {
			
            if (dataX[i]>maxx) {
                dataY[i] = (2* maxx) - dataX[i];
            } else {
                dataY[i] = dataX[i];
            }
        }
		
        for (int i=0;i<dataY.length;i++) {
			
            if (dataY[i]<minx) {
                dataY[i] = (2* minx) - dataY[i];
            } else {
                dataY[i] = dataX[i];
            }
			
            if (dataX[i]>maxx) {
                dataY[i] = (2* maxx) - dataY[i];
            } else {
                dataY[i] = dataY[i];
            }
			
        }
        return dataY;
		
    }
	
	
    public static double[] convolve(double[] data, double[] operator, double[] output){
		
        int dataLen = data.length;
        int operatorLen = operator.length;
		
        for(int i = 0;i < dataLen-operatorLen+1;i++){
            output[i] = 0;
			for(int j = operatorLen-1;j >= 0;j--){
				output[i] += data[i+j]*operator[j];
			}
        }
        return output;
    }
	
    public static double[] extend(double[] data) {
		
        double[] newdata = new double[data.length+2];
        
        newdata[0] = data[0];
        for (int i=0;i<data.length;i++) {
            newdata[i+1] = data[i];
        }
		
        newdata[newdata.length-1] = data[data.length-1];
		
        return newdata;
    }
	
	
    public static double[] Interpol(double[] X, double med) {
		
        double[] XX = Absol(X);
        double quant = quantile(XX, 0.68);
        
        for (int i=0;i<XX.length;i++) {          
            if (XX[i] > quant) {
                XX[i] = med;
				
            } else {
                XX[i] = X[i];
				
            }          
        }
        
        return XX;     
    }
	
	public static Object[][] pWave(Object[][] data) {
		
		
 		Object[] chrs = new Object[data.length];
 		for(int i=0;i<data.length;i++) {
 			chrs[i] = data[i][1]; 
 		}
 		
 		ArrayList uni = new ArrayList();
 		java.util.List lis = Arrays.asList(chrs);
		for (Object x : lis) {
    		if (!uni.contains(x)) {
        		uni.add(x);
    		}
   		}
   		
   		ArrayList[] ratio = new ArrayList[uni.size()];
   		ArrayList[] calls = new ArrayList[uni.size()];
   		ArrayList[] starts = new ArrayList[uni.size()];
   		ArrayList[] stops = new ArrayList[uni.size()];
   		
   		//TRY:	index out of bounds
   		/*int ma = 0;
   		for(int i=0;i<uni.size();i++) {
   			int th = Integer.parseInt(""+uni.get(i));
   			if(th>ma) {
   				ma=th;
   			}
   		}
   		
   		ArrayList[] ratio = new ArrayList[ma];
   		ArrayList[] calls = new ArrayList[ma];
   		ArrayList[] starts = new ArrayList[ma];
   		ArrayList[] stops = new ArrayList[ma];*/
   		
   		for(int i=0;i<ratio.length;i++) {
   			calls[i] = new ArrayList();
   			ratio[i] = new ArrayList();
   			starts[i] = new ArrayList();
   			stops[i] = new ArrayList();
   		}
   		for(int i=0;i<data.length;i++) {
   			int loc = Integer.parseInt(""+data[i][1]) - 1;
   			calls[loc].add(data[i][0]);
   			ratio[loc].add(data[i][2]);
   			starts[loc].add(data[i][3]);
   			stops[loc].add(data[i][4]);
   		}
   		
   		int lpin = 0;
		
   		for(int i=0;i<uni.size();i++) {
   			
   			Object[] rats = new Object[ratio[i].size()];
   			rats = ratio[i].toArray(rats);
   			double m = mean(rats);
   			
   			Object[] ca = new Object[calls[i].size()];
   			ca = calls[i].toArray(ca);
   			Object[] sta = new Object[starts[i].size()];
   			sta = starts[i].toArray(sta); 
   			Object[] sto = new Object[stops[i].size()];
   			sto = stops[i].toArray(sto); 
   			
   			ca = runWave(i, rats,ca, sta, sto);
			
			for(int k =0; k< ratio[i].size(); k++) {
				ratio[i].set(k, ca[k]);
				data[lpin][2] = ca[k];
				lpin++;
			}
   			
   		}
		
   		
   		return(data);
    }
	
	
	public static Object[] runWave(int chr, Object[] ratios, Object[] calls, Object[] sta, Object[] sto) {
		
    	int pin = 0;
    	int pin2 = 0;
    	int index = 0;
		int check = 0;
		
		System.out.println(chr + "\t" + calls.length + "\t" + ratios.length);
    	boolean seen = false;
    	ArrayList meds = new ArrayList();
    	double m1 = mean(ratios);
    	
		Object[] refCalls = new Object[calls.length];
		for(int i=0;i<refCalls.length;i++) {
			refCalls[i] = m1;
		}
		
    	for(int i=0;i<ratios.length;i++) {
  			
  			index = i;
			seen = true;
			
    		if(Integer.parseInt("" + calls[i]) == 1) {    			
    			
    			if (seen & check == 0) {
    				pin2 = i;
    			}
				check++;
    			pin = 1;
    			while(Integer.parseInt("" + calls[index]) == 1 & index <calls.length) {
    				index++;
    				//System.out.println(calls.length + "\t" + index);
    				if (index>=calls.length) {
    					pin = 0;
    					seen = false;
    					check=0;
    					break;
    				}
    			}
    			
    		} else {
    			pin = 0;
    			seen = false;
    			check=0;
    		}
    		
   			if (index>=calls.length) { break; }
    		if (pin == 1 & i==index-1 & i != 0 & pin2 != i) { 
				//System.out.println(i + "\t" + pin2);
    		 	Object[] temp = new Object[i-pin2-1];
   				System.arraycopy(ratios, pin2, temp, 0, temp.length);
   				double m2 = mean(temp);
				for(int k = pin2;k<temp.length;k++) {
					refCalls[k] = m2;
				}
				
    		}

    	}
		
		double[] ddd = new double[ratios.length];
    	
    	for(int i=0;i<refCalls.length;i++) {
    			ddd[i] = Double.parseDouble(""+ ratios[i]) - Double.parseDouble(""+refCalls[i]);
    	}
		
		double[] D = Interpol(ddd, 0);
		double[] my = dtWaveXfm(D, 1);
		Object[] R = new Object[ratios.length];		
		
		for (int i=0;i<ddd.length;i++) {
         	R[i] = Double.parseDouble(""+ratios[i]) -my[i];
		}

    	return R;
    }
	
    
    public static Object[][] loopRep(ArrayList[] report, ArrayList[] data) {
    	
    	int dchr = 0;
    	int dstart = 1;
    	int dstop = 1;
    	int pin1 = 0;
    	int pin2 = 0;
    	
    	Object[][] isCall = new Object[data[0].size()][5];
    	for(int i=0;i<isCall.length;i++) {
    		isCall[i][0] = 0;
    		isCall[i][1] = Integer.parseInt(""+ data[0].get(i));
    		isCall[i][2] = Double.parseDouble(""+ data[3].get(i));
    		isCall[i][3] = Integer.parseInt(""+ data[1].get(i));
    		isCall[i][4] = Integer.parseInt(""+ data[2].get(i));
    	}
    	
    	for(int i=0;i<report[0].size();i++) {
			
 			boolean go = true;
    		int rchr = Integer.parseInt(""+ report[0].get(i));
    		int rstart = Integer.parseInt(""+ report[1].get(i));
    		int rstop = Integer.parseInt(""+ report[2].get(i));
			
    		while(go) {
				
    			dchr = Integer.parseInt(""+ data[0].get(pin1));
    			dstart = Integer.parseInt(""+ data[1].get(pin1));
				
    			if (rchr == dchr & rstart == dstart) { 
    				
    				boolean go2 = true;
    				pin2 = pin1;
    				
					while(go2) {				
						
						//TRY: change due to index error!
						if (pin2 >= data[0].size()-1) { pin2 = isCall.length-1; go2 = false; break; }
						//if (pin2 >= data[0].size()) { pin2 = isCall.length-1; go2 = false; break; }
						dstop = Integer.parseInt(""+ data[2].get(pin2));	
						
						if (rchr == dchr & rstop == dstop) {
							go2 = false;
						}
						
						pin2++;
					}
					
    				go=false;
    			}

    			pin1++;
				
    		}
    		
 			for(int k=pin1;k<=pin2;k++) {
 				isCall[k][0] = 1;
 			}
		}
    	
    	return(isCall);    	
    }
	
	public static ArrayList[] getData(String filename) {
		
     	ArrayList[] data = new ArrayList[4];
		for (int i=0;i<data.length;i++) {
			data[i] = new ArrayList();
		}
		
		filename.trim();
		String inLine;
		BufferedReader infile = null;
		
		try {
			infile = new BufferedReader(new FileReader (filename));
			
			while ((inLine=infile.readLine()) != null) {
				
				String[] Ldata = inLine.split("\t");
				data[0].add(Ldata[0]);
				data[1].add(Ldata[1]);
				data[2].add(Ldata[2]);
				data[3].add(Ldata[3]);
				
			}
		}
		catch (FileNotFoundException ex) {
			System.out.println("File not found: " + filename);
		}
		catch (IOException ex) {
			System.out.println(ex.getMessage());
		}
		finally {
			try {
				if (infile != null) infile.close();
			}
			catch (IOException ex) {
				System.out.println(ex.getMessage());
			}
		}
		
		return data;
		
    }
    
    public static ArrayList[] getReport(String filename) {
		
     	ArrayList[] data = new ArrayList[4];
		for (int i=0;i<data.length;i++) {
			data[i] = new ArrayList();
		}
		
		filename.trim();
		String inLine;
		BufferedReader infile = null;
		
		try {
			infile = new BufferedReader(new FileReader (filename));
			
			while ((inLine=infile.readLine()) != null) {
				
				String[] Ldata = inLine.split("\t");
				data[0].add(Ldata[0]);
				data[1].add(Ldata[1]);
				data[2].add(Ldata[2]);
				data[3].add(Ldata[3]);
				
			}
		}
		catch (FileNotFoundException ex) {
			System.out.println("File not found: " + filename);
		}
		catch (IOException ex) {
			System.out.println(ex.getMessage());
		}
		finally {
			try {
				if (infile != null) infile.close();
			}
			catch (IOException ex) {
				System.out.println(ex.getMessage());
			}
		}
		
		return data;
		
    }
	
	public static void OutFILE(Object[][] calls, String filename) {
		
		PrintWriter pw = null;
		File output = new File(filename);    
		
		try {
			pw = new PrintWriter(new FileOutputStream(output), true);
			
			for(int i=0;i<calls.length;i++) {
				pw.print(calls[i][1] + "\t" + calls[i][3] + "\t" + calls[i][4] + "\t" + calls[i][2] + "\n");
			}
			
		}
		catch (IOException ex) {
			System.out.println(ex.getMessage());
		}
		finally {
			if (pw != null) pw.close();
		}
	}
	
	
    
   public static void main(String[] args) {
	
		String file ="";
		String report ="";
		
		for (int u=0;u<args.length; u++) {
			
			int pin = u+1;
            String s = args[u]; 
			
			if (s.equals("-f")) { 
				file = args[pin];
			} 
			if (s.equals("-r")) { 
				report = args[pin];
			} 
		}
		
		ArrayList[] data = getData(file);
		ArrayList[] rep = getReport(report);
		Object[][] calls = loopRep(rep, data);
		System.out.println("Called");
		calls = pWave(calls);
	    OutFILE(calls, file);
	    Clean.ExeShellCommand("rm " + report);
	}
    

}






