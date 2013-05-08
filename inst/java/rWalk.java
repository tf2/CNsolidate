import java.lang.Math.*;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;

public class rWalk {

	public static int ArraySizer(String filename) {
     filename.trim(); String inLine; BufferedReader infile = null; int pin =0;
     try {
         infile = new BufferedReader(new FileReader (filename));
         while ((inLine=infile.readLine()) != null) {
            pin++;
         }
     } catch (FileNotFoundException ex) { System.out.println("File not found: " + filename); }
     catch (IOException ex) { System.out.println(ex.getMessage()); }
     finally { try { if (infile != null) infile.close(); }
         catch (IOException ex) { System.out.println(ex.getMessage());  }
     }
    return pin;
}

    public static double[] In(String filename, int size) {
    	filename.trim(); String inLine; BufferedReader infile = null;
     	double [] elements  = new double [size]; int pin =0;
     	try {
         	infile = new BufferedReader(new FileReader (filename));
         	while ((inLine=infile.readLine()) != null && pin <size) {
            	String[] seq = inLine.split("\t");
            	String el = seq[0];   
            	if( el.equals("NA") || el.equals("NaN") ) el = "" +0;
            	elements[pin] = Double.parseDouble(el);
            	pin++;
         	}
     	} catch (FileNotFoundException ex) { System.out.println("File not found: " + filename); }
     	catch (IOException ex) { System.out.println(ex.getMessage()); }
     	finally { try { if (infile != null) infile.close(); }
         catch (IOException ex) { System.out.println(ex.getMessage()); }
     	}
    return elements;
	}

	public static double var(double[] data) {
        int n = data.length; int d = n-1;
        double[] xbar = new double[data.length]; 
        double m = mean(data); double var = 0;
        for (int i=0;i<n;i++) {
            double deviation = data[i] - m;
            var += Math.pow(deviation,2);
        }
        return var/d;
    }

	public static double mean(double[] m) {
      double men = 0; int l = m.length;  
      for (int i = 0 ; i<m.length ; i++) men = men + m[i];
      return men / l;
   	}

	public static double median(double[] values) {    
      double[] copy = new double[values.length]; System.arraycopy(values, 0, copy, 0, copy.length);
      java.util.Arrays.sort(copy); int middle = copy.length/2;
      if (copy.length%2 == 1) return copy[middle];
	  else return (copy[middle-1] + copy[middle]) / 2.0;
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
  	
  	
  	public static void rWalking(double[] values) {
  		int var_win = 3; double thred_set = 0.05; double thred_bound = 0.1; int thred_pt_num = 2; 	
  		double[] data_filter = runMedian(values, 3); double[] data_filter_int = new double[data_filter.length]; double men = mean(data_filter);
    	for(int i=0;i<data_filter.length;i++) { data_filter[i]=data_filter[i]-men;  data_filter_int[i] = data_filter[i]; }
   		for (int i=1;i<data_filter.length;i++)  data_filter_int[i] = data_filter_int[i-1] + data_filter[i];	
		int DL = data_filter_int.length; int win_len = var_win; int window_overlap_size = 1; int S = (int)Math.floor((DL-win_len)/window_overlap_size+1)+1;
		double[] keep=new double[S];
		for (int j=1;j<S-1;j++) {	
			int st = (j - 1)*window_overlap_size; int so = (j - 1)*window_overlap_size + win_len;	
			double[] copy = new double[so-st]; System.arraycopy(data_filter_int, st, copy, 0, copy.length);
			keep[j+(win_len-1)/2] = var(copy);

		} double[] keep_new=new double[keep.length];
		for (int i = (win_len + 1)/2; i< keep.length-(win_len - 1)/2; i++) {
			double[] va = new double[2]; System.arraycopy(data_filter, i-1, va, 0, va.length); double v1 = var(va);
    		va = new double[2]; System.arraycopy(data_filter, i, va, 0, va.length); double v2 = var(va);	
    		if (v1 >= v2) { keep_new[i] = data_filter[i+1]; } else { keep_new[i] = data_filter[i-1]; }
		}
		double[] result=new double[keep_new.length];
		for (int i=0;i<keep.length;i++) {
			result[i]=keep_new[i]*keep[i];
		} double[] result_new=new double[result.length];
		for (int i=0;i<keep.length;i++) {
			result_new[i]=result[i];
		}
		for (int i = 1; i<result.length;i++) { if (Math.abs(result[i]) > thred_set)  result_new[i] = keep_new[i]; }

		double[] result_change=new double[result_new.length];
		for (int i=0;i<keep.length;i++) result_change[i]=result_new[i];

		boolean seen = false; int pin=0; double t1=0;
		for (int i = 0; i< result_change.length-1;i++) {
			if(result_change[i] > thred_bound) {
    			t1 = result_change[i]; pin=i;
      			while(result_change[pin]>thred_bound) {
      				pin++;
      			} seen = true;
     		}
   	 		if(seen) {
     			for(int j=i;j<pin;j++) {
        			result_change[j] = (t1 + result_change[pin-1])/2;
        		}
    		}
		}
		for (int i=0;i<result.length;i++) {
			System.out.println(result[i] + "\t" + result_new[i] + "\t" + result_change[i]); 
		}
  	}
  	
  	public static void main(String[] args) {
  		//double[] values =  {-0.2807,-0.5492, 0.2946, -0.1036, 0.1202, 0.0956, 0.0393, -0.2321, -0.2170, -0.0655, -0.2802, -0.2458, 0.4233, 0.0484, -0.2885, 0.0227, -0.2650, -0.3235,-0.0355,0.2859};   
  		double[] values = In("example_data.txt", ArraySizer("example_data.txt"));
  		rWalking(values);

	
    }
    
  	
}