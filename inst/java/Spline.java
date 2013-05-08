/*
Author:	Tomas William Fitzgerald
All rights reserved.
*/

import java.util.regex.*;
import java.io.*;
import java.io.FileNotFoundException;
import java.util.Arrays;
import java.util.Vector;
import java.util.logging.Level;
import java.util.logging.Logger;


public class Spline {

     public static String NoiseT; public static double thes = 0.68;
     public static double fact = 4.5; public static int knots = 1000; public static int iter = 5;

    public Spline() {

    }
    
    public static double quantile(double[] values, double quantile) {
        if(values == null) throw new IllegalArgumentException("Values cannot be null.");
        if(quantile < 0.0 || quantile > 1.0) throw new IllegalArgumentException("Quantile must be between 0.0 and 1.0");
        double[] copy = new double[values.length]; System.arraycopy(values, 0, copy, 0, copy.length); java.util.Arrays.sort(copy); int index = (int) (copy.length * quantile);
    return copy[index];
    }
    
    public static double median(double[] values) {
		double[] copy = new double[values.length]; System.arraycopy(values, 0, copy, 0, copy.length);
      	java.util.Arrays.sort(copy); int middle = copy.length/2;
      	if (copy.length%2 == 1) { return copy[middle]; } else { return (copy[middle-1] + copy[middle]) / 2.0; }
	}
    
     public static int Sizer(String filename) {
     	Pattern p = Pattern.compile("chr\\w\\w?:\\d+-\\d");
     	filename.trim(); String inLine;  String outline = "";  BufferedReader infile = null; int len=0;      
     	try { infile = new BufferedReader(new FileReader (filename));
        	 while ((inLine=infile.readLine()) != null) {  Matcher m = p.matcher(inLine); if (m.find()) len++; }
     	}
     	catch (FileNotFoundException ex) { System.out.println("File not found: " + filename); }
     	catch (IOException ex) { System.out.println(ex.getMessage()); }
     	finally {  try { if (infile != null) infile.close(); } catch (IOException ex) { System.out.println(ex.getMessage()); } }    
    return len;       
    }
    
     public static int Sizer2(String filename) {
     	filename.trim(); String inLine; BufferedReader infile = null; int len=0;      
     	try { infile = new BufferedReader(new FileReader (filename));
        	 while ((inLine=infile.readLine()) != null) { len++; }
     	}
     	catch (FileNotFoundException ex) { System.out.println("File not found: " + filename); }
     	catch (IOException ex) { System.out.println(ex.getMessage()); }
     	finally {  try { if (infile != null) infile.close(); } catch (IOException ex) { System.out.println(ex.getMessage()); } }    
    return len;       
    }
    
    public static double[][] Input(String file, int size) {      
     	DataInputStream dis = null; String record = null; 
     	int recCount = 0; int probe_int = 0; int red_int =0; int green_int =0;  int loc_int =0; int ratio_int =0;
        Pattern p = Pattern.compile("chr\\w\\w?:\\d+-\\d"); String da = ""; int pin =0; double[][] ret_vec = new double[size][2];     
        try { 
 			File f = new File(file); FileInputStream fis = new FileInputStream(f); 
           	BufferedInputStream bis = new BufferedInputStream(fis); dis = new DataInputStream(bis);  
           		while ( (record=dis.readLine()) != null ) { 
              		recCount++; 
             		if (recCount == 10){ String[] seq = record.split("\t");
              			for (int i =0; i<seq.length; i++ ) {
              				String temp = seq[i];
                 			if (temp.equals("ProbeName")) probe_int = i;
                			else if (temp.equals("LogRatio")) ratio_int = i;
                			else if (temp.equals("rProcessedSignal")) red_int = i;
                			else if (temp.equals("gProcessedSignal")) green_int = i;
                 			//else if (temp.equals("rBGSubSignal")) red_int = i;
                 			//else if (temp.equals("gBGSubSignal")) green_int = i;
                 			else if (temp.equals("SystematicName")) loc_int = i;                  
              			}
           			}
              		if (recCount > 10) {
                  		String[] data = record.split("\t"); Matcher matc= p.matcher(record); 
                  
                  		if (matc.find()) {
                      		ret_vec[pin][0] = Double.parseDouble(data[red_int]); ret_vec[pin][1] = Double.parseDouble(data[green_int]);
                  			pin++;
                		}                 
             	 	}
           		}
       			} catch (IOException e) { System.out.println("Uh oh, got an IOException error!" + e.getMessage()); 
        		} finally { if (dis != null) { try { dis.close(); } catch (IOException ioe) {}
           		} 
        	}
      return ret_vec;
      }

 	public static void Output(String file, String filename, double[] red, double[] green) {
    PrintWriter pw = null; File output = new File(filename); double[] flag = new double[red.length];
    if(output.exists()) System.out.println("File already exists - overwriting!!!");
     	DataInputStream dis = null; String record = null;
        int recCount = 0; int loc_int =0; int pro_int =0; int pin =0; int raterr_int=0;
        Pattern p = Pattern.compile("chr\\w\\w?:\\d+-\\d");
        try {
           File f = new File(file); FileInputStream fis = new FileInputStream(f); BufferedInputStream bis = new BufferedInputStream(fis);
           dis = new DataInputStream(bis); pw = new PrintWriter(new FileOutputStream(output), true);
           while ( (record=dis.readLine()) != null ) {
           recCount++;
        	if (recCount == 10){
            	String[] seq = record.split("\t");
              for (int i =0; i<seq.length; i++ ) {
              	String temp = seq[i];
                if (temp.equals("SystematicName")) loc_int = i;
                else if (temp.equals("ProbeName")) pro_int = i;
                else if (temp.equals("LogRatioError")) raterr_int = i;
              }
           }
              else if (recCount > 10) {
                  String[] data = record.split("\t");
                  Matcher matc = p.matcher(record);
                  if (matc.find()) {
                  	StringBuilder sb = new StringBuilder(512); boolean first = true; String SyName = data[loc_int];
                    String[] SyArray = SyName.split(":"); String SyTemp = SyArray[1]; String[] SyLocat = SyTemp.split("-");                   	
                      	int start = Integer.parseInt(SyLocat[0]); int end = Integer.parseInt(SyLocat[1]); int temp_pos;
                      	if (start > end) { temp_pos = end; end = start; start = temp_pos; }
                      	if (red[pin]+green[pin]<100) flag[pin] = 1;
                      	else flag[pin] = 0;
                      	double LogRat2 = Math.log(red[pin]/ green[pin]) / Math.log(2);
                      	double err = 1-Double.parseDouble(""+data[raterr_int]);
                       	sb.append(SyArray[0] + "\t" + start + "\t" + end + "\t" + LogRat2 + "\t" + err + "\t" + red[pin] + "\t" + green[pin] + "\t" + data[pro_int] + "\n");
                       	pw.print(sb); pin++;
                }
              }
           }
        	} catch (IOException e) { System.out.println("Uh oh, got an IOException error!" + e.getMessage());
        	} finally { if (dis != null) { try { pw.close(); dis.close();
	      } catch (IOException ioe) {}
          }
      	}
     }

	public static void Output2(String file, String filename, double[] raw_red, double[] raw_green, double[] red, double[] green) {
    PrintWriter pw = null; File output = new File(filename); double[] flag = new double[red.length];
    if(output.exists()) System.out.println("File already exists - overwriting!!!");
     	DataInputStream dis = null; String record = null;
        int recCount = 0; int loc_int =0; int pro_int =0; int pin =0; int raterr_int=0;
        Pattern p = Pattern.compile("chr\\w\\w?:\\d+-\\d");
        try {
           File f = new File(file); FileInputStream fis = new FileInputStream(f); BufferedInputStream bis = new BufferedInputStream(fis);
           dis = new DataInputStream(bis); pw = new PrintWriter(new FileOutputStream(output), true);
           while ( (record=dis.readLine()) != null ) {
           recCount++;
        	if (recCount == 10){
            	String[] seq = record.split("\t");
              for (int i =0; i<seq.length; i++ ) {
              	String temp = seq[i];
                if (temp.equals("SystematicName")) loc_int = i;
                else if (temp.equals("ProbeName")) pro_int = i;
                else if (temp.equals("LogRatioError")) raterr_int = i;
              }
           }
              else if (recCount > 10) {
                  String[] data = record.split("\t");
                  Matcher matc = p.matcher(record);
                  if (matc.find()) {
                  	StringBuilder sb = new StringBuilder(512); boolean first = true; String SyName = data[loc_int];
                    String[] SyArray = SyName.split(":"); String SyTemp = SyArray[1]; String[] SyLocat = SyTemp.split("-");                   	
                      	int start = Integer.parseInt(SyLocat[0]); int end = Integer.parseInt(SyLocat[1]); int temp_pos;
                      	if (start > end) { temp_pos = end; end = start; start = temp_pos; }
                      	if (red[pin]+green[pin]<100) flag[pin] = 1;
                      	else flag[pin] = 0;
                      	double LogRat2 = Math.log(red[pin]/ green[pin]) / Math.log(2);
                      	double err = 1-Double.parseDouble(""+data[raterr_int]);
                       	sb.append(SyArray[0] + "\t" + start + "\t" + end + "\t" + LogRat2 + "\t" + err + "\t" + raw_red[pin] + "\t" + raw_green[pin] + "\t" + red[pin] + "\t" + green[pin] + "\t" + data[pro_int] + "\n");
                       	pw.print(sb); pin++;
                }
              }
           }
        	} catch (IOException e) { System.out.println("Uh oh, got an IOException error!" + e.getMessage());
        	} finally { if (dis != null) { try { pw.close(); dis.close();
	      } catch (IOException ioe) {}
          }
      	}
     }

	public static double[][] InputBED(String file, int size) {      
     	DataInputStream dis = null; String record = null; 
		int pin =0; double[][] ret_vec = new double[size][2];     
        try { 
 			File f = new File(file); FileInputStream fis = new FileInputStream(f); 
           	BufferedInputStream bis = new BufferedInputStream(fis); dis = new DataInputStream(bis);  
			while ( (record=dis.readLine()) != null ) {
				String[] data = record.split("\t");
				ret_vec[pin][0] = Double.parseDouble(data[5]); ret_vec[pin][1] = Double.parseDouble(data[6]);
				//System.out.println(data[5] + "\t" + data[6]);
				pin++;
			}
		} catch (IOException e) { System.out.println("Uh oh, got an IOException error!" + e.getMessage()); 
		} finally { if (dis != null) { try { dis.close(); } catch (IOException ioe) {}
		} 
		}
		return ret_vec;
	}
	
	public static void OutputBED(String file, String filename, double[] red, double[] green) {
		PrintWriter pw = null; File output = new File(filename); double[] flag = new double[red.length];
		if(output.exists()) System.out.println("File already exists - overwriting!!!");
     	DataInputStream dis = null; String record = null; int pin=0;
        try {
			File f = new File(file); FileInputStream fis = new FileInputStream(f); BufferedInputStream bis = new BufferedInputStream(fis);
			dis = new DataInputStream(bis); pw = new PrintWriter(new FileOutputStream(output), true);
			while ( (record=dis.readLine()) != null ) {
					String[] data = record.split("\t");
						StringBuilder sb = new StringBuilder(512); boolean first = true;
                      	if (red[pin]+green[pin]<100) flag[pin] = 1;
                      	else flag[pin] = 0;
                      	double LogRat2 = Math.log(red[pin]/ green[pin]) / Math.log(2);
                       	sb.append(data[0] + "\t" + data[1] + "\t" + data[2] + "\t" + LogRat2 + "\t" + data[4] + "\t" + red[pin] + "\t" + green[pin] + "\t" + flag[pin] + "\n");
                       	pw.print(sb); pin++;
			}
		} catch (IOException e) { System.out.println("Uh oh, got an IOException error!" + e.getMessage());
		} finally { if (dis != null) { try { pw.close(); dis.close();
		} catch (IOException ioe) {}
		}
      	}
	}
	
	public static void OutputBED2(String file, String filename, double[] raw_red, double[] raw_green, double[] red, double[] green) {
		PrintWriter pw = null; File output = new File(filename); double[] flag = new double[red.length];
		if(output.exists()) System.out.println("File already exists - overwriting!!!");
     	DataInputStream dis = null; String record = null; int pin=0;
        try {
			File f = new File(file); FileInputStream fis = new FileInputStream(f); BufferedInputStream bis = new BufferedInputStream(fis);
			dis = new DataInputStream(bis); pw = new PrintWriter(new FileOutputStream(output), true);
			while ( (record=dis.readLine()) != null ) {
					String[] data = record.split("\t");
						StringBuilder sb = new StringBuilder(512); boolean first = true;
                      	if (red[pin]+green[pin]<100) flag[pin] = 1;
                      	else flag[pin] = 0;
                      	double LogRat2 = Math.log(red[pin]/ green[pin]) / Math.log(2);
                       	sb.append(data[0] + "\t" + data[1] + "\t" + data[2] + "\t" + LogRat2 + "\t" + data[4] + "\t" + raw_red[pin] + "\t" + raw_green[pin] + "\t" + red[pin] + "\t" + green[pin] + "\t" + flag[pin] + "\n");
                       	pw.print(sb); pin++;
			}
		} catch (IOException e) { System.out.println("Uh oh, got an IOException error!" + e.getMessage());
		} finally { if (dis != null) { try { pw.close(); dis.close();
		} catch (IOException ioe) {}
		}
      	}
	}
	
	
    public static int Isearch(double[] data, double point) {
        int ipin = 0; int i = 0;  int pin = 0;   
        while(ipin == 0 && i < data.length-1) {
        	double t = Math.abs(data[i] - point); double tt = Math.abs(data[i+1]-point);
            if (t < tt) { ipin = 1; pin = i;  }
            if (i == data.length-2) { ipin = 1; pin = i-1; }
            i++;
        }
    return(pin);
    }

    public static double[][] FE(String file) {
        double[][] my = Input(file, Sizer(file));
   	return (my);
    }
	
	public static double[][] BED(String file) {
        double[][] my = InputBED(file, Sizer2(file));
		return (my);
    }

	public static double[] RED(double[] x, double[] y, int knots) {
        double[] result =eSpline(x, y, iter, knots);
        x = null; y = null;
    return (result);
    }

    public static double[] eSpline(double[]x, double[] y, int iterate, int knots) {
       	double[] result = null; result = new double[x.length]; int nu=x.length; double e=0; int colum=2; int len=x.length;
       	double number = Math.round(len/knots); number = Math.max(number, 100); double[] pin = new double[(int)number];
     	for (int i=0;i<pin.length;i++) pin[i] = (i+1) / number;
       	int[] ind = new int[(int)number];
        for (int i=0;i<ind.length;i++) ind[i] = (int) (pin[i] * len);
       	int Off = Math.round(ind[0]/iterate); int[] Offa = new int[iterate]; Offa[0] = 0;
       	for (int i=1;i<iterate;i++) Offa[i] = Off;
		double[] xx = new double[x.length];
			for (int i=0 ; i<x.length;i++) {
				if (x[i] < 0) xx[i] = Math.abs(xx[i]);
				else xx[i] = x[i];
           }
       	double[] sY = new double[x.length];
         	for (int u = 0; u<sY.length; u++) {
            	if(x[u]<0) x[u] = Math.abs(x[u]);
             	if(y[u]<0) y[u] = Math.abs(y[u]);
             	sY[u] = Math.exp((Math.log(x[u]) + Math.log(y[u])) / 2);
			 	if (sY[u] < 0) sY[u] = 1;
         	}
       	double[] sX = new double[x.length];
       	for (int i=0;i<sX.length;i++) sX[i] = xx[i];
       	Arrays.sort(sX); Arrays.sort(sY);
       	
       	//start offset!!
       	for (int p = 0; p<Offa.length;p++) {
           	for (int i=0;i<ind.length;i++) ind[i] = ind[i] - Offa[p];
       	   	double[] xPol = new double[(int)number];
       	   	for(int i=0;i<xPol.length;i++) xPol[i] = sX[(int)ind[i]-1];
		   	xPol[xPol.length-1] = sX[len-1];
       		double[] yPol = new double[(int)number];
	   		for (int i=0;i<yPol.length;i++) yPol[i] = sY[(int)ind[i]-1];
			yPol[yPol.length-1] = sY[len-1];
         	for(int i = 0; i < xPol.length; i++) {
            	for(int j = i+1; j < xPol.length; j++) {
                	if(xPol[j] < xPol[i]) { double temp = xPol[i]; xPol[i] = xPol[j]; xPol[j] = temp; } 
                	if(yPol[j] < yPol[i]) { double temp2 = yPol[i]; yPol[i] = yPol[j]; yPol[j] = temp2; }
            	}
         	}
        	Vector<Double> vector1 = new Vector(); Vector<Double> vector2 = new Vector();
     		int count=0; int pin2=0; double tempA=0; double tempB=0;
     			for (int i=0;i<xPol.length;i++) {
           		double tempX =xPol[i]; double tempY =yPol[i];
            		if (i == xPol.length-1) {
                  		double ave1 = 0; double ave2=0; tempX = xPol[i]; tempY =yPol[i];       	
                		if (count == 0) {
                 			count =1; ave1 = tempX; ave2 = tempY;
                		} else if (count > 0) {
                 			ave1 = tempA / count; ave2 = tempB / count;
                		} vector2.add(pin2, ave1); vector1.add(pin2, ave2); tempA =0; count=0;
                		} else if (i < xPol.length-1) {
                    		int k = i+1;    
                    		if (xPol[i] == xPol[k]) { 
                    			tempA =  tempA + tempX; tempB = tempB + tempY; count++;
                    		} else {
                    			double ave1=0; double ave2=0;
                    		if (count == 0) { 
                    			ave1 = tempX; ave2 = tempY;
                    		} else if(count>0) { 
                    			ave1 = tempA / count; ave2 = tempB / count;
                    		}
               			vector2.add(pin2, ave1); vector1.add(pin2, ave2);
               			tempA =0; tempB =0; count=0; pin2++;
               			}
               		}
           		}
         	xPol = null; yPol = null; System.gc();
         	int size = vector1.size(); double[] uniY = new double[size]; double[] uniX = new double[size];
         	for (int i=0;i<uniX.length;i++) {
             	double temp = vector1.get(i); uniY[i] = temp; temp = vector2.get(i); uniX[i] = temp;
         	} vector1 = null; vector2 = null; System.gc();
     		
     		// SPLINE
       		double[] xxx = JsplineFun(uniX,uniY,xx);
       		for (int i =0;i<result.length;i++) result[i] += xxx[i] /iterate;
	 		xxx=null; uniX = null; uniY = null; System.gc();
     
     	} // end offset!!	   
	   	xx=null; sX = null; sY = null; System.gc();
	return(result);
	}

   
    public static double[] eINTER(double[] red, double[] green) {      
    	int size = red.length; double[] lo2 = new double[size];
        for (int i=0;i<red.length;i++) {
            if (red[i] <= 0) red[i] = Math.abs(red[i]);
            if (green[i] <= 0) green[i] = Math.abs(green[i]);
        }  
        for (int i=0; i<size;i++) lo2[i] = Math.log(red[i]/green[i]) / Math.log(2);
        double med =median(lo2);
        for (int i=0; i<size;i++) lo2[i] = Math.abs(lo2[i] - med);  
        double rp=0; rp=quantile(lo2, thes);
        rp = rp * fact;
        Vector<Double> FvetR = new Vector(); Vector<Double> FvetG = new Vector();  Vector<Double> IvetR = new Vector();
        int pin=0;
        for (int i=0;i<size;i++) {
            if (lo2[i] <= rp) { FvetR.add(pin, red[i]); FvetG.add(pin, green[i]); pin++; }
        }
        pin=0;
        for (int i=0;i<size;i++) { if (lo2[i] > rp) { IvetR.add(pin, red[i]); pin++; } }
        red = null; green=null; System.gc();
        int size1 = FvetR.size(); int size2 = IvetR.size();
        double[] fX = new double[size1]; double[] fY = new double[size1];
        for (int i=0;i<fX.length;i++) { double temp = FvetR.get(i); fX[i] = temp; temp = FvetG.get(i); fY[i] = temp; }
        FvetR = null; FvetG = null; System.gc();
        double[] fit = RED(fX,fY, knots);
        fY=null; System.gc();
        double[] iX = new double[size2];
        for (int i=0;i<iX.length;i++) { double temp = IvetR.get(i); iX[i] = temp; }
        double[] ifit = INTER(iX, fX, fit);
        iX = null; fX = null; System.gc();
        double[] result = new double[size];
        int p1 = 0; int p2 = 0;
        for (int i=0 ; i<size;i++) {
            if (lo2[i] <= rp) { result[i] = fit[p1]; p1++; } else if (lo2[i] > rp) { result[i] = ifit[p2]; p2++; }
        }
        fit = null; ifit = null; System.gc();
    return (result);
    }

    public static double[] JsplineFun(double[] xX, double[] yY, double[] testX) {
    	int sN = xX.length; int sNN = testX.length; int i, j, k, l;  
    	double t, ul, dx, tmp; int nm1 = sN - 1; int n_1 = sNN - 1;
    	double[] sX = new double[sN+1]; double[] sY = new double[sN+1]; double[] sB = new double[sN+1];
    	double[] sC = new double[sN+1]; double[] sD = new double[sN+1]; double[] sV = new double[sNN];
    	double[] retArr = new double[sNN];

    	for (int p =0;p<sNN;p++) { retArr[p] = testX[p]; sV[p] = testX[p]; }

    	for (int p=sN;p>0;p--) {
        	sX[p] = xX[p-1]; sY[p] = yY[p-1];
        	sB[p] = 0; sC[p] = 0; sD[p] = 0;
    	} sX[0] = 0; sY[0] = 0; sB[0] = 0; sC[0] = 0; sD[0] = 0;
    
    	Arrays.sort(sX); Arrays.sort(sY);
    	sD[1] = sX[2] - sX[1]; sC[2] = (sY[2] - sY[1])/sD[1];
    	for( i=2 ; i<sN ; i++) {
			sD[i] = sX[i+1] - sX[i];
			sB[i] = 2.0 * (sD[i-1] + sD[i]);
			sC[i+1] = (sY[i+1] - sY[i])/sD[i];
			sC[i] = sC[i+1] - sC[i];
    	}

    	for(i=3 ; i<sN ; i++) {
			t = sD[i-1]/sB[i-1];
			sB[i] = sB[i] - t*sD[i-1];
			sC[i] = sC[i] - t*sC[i-1];
    	} sC[nm1] = sC[nm1]/sB[nm1];
    
    	for(i=sN-2 ; i>1 ; i--) {
			sC[i] = (sC[i]-sD[i]*sC[i+1])/sB[i];
    	} sC[1] = sC[sN] = 0.0; sB[sN] = (sY[sN] - sY[sN-1])/sD[sN-1] + sD[sN-1]*(sC[sN-1]+ 2.0*sC[sN]);

    	for(i=1 ; i<=nm1 ; i++) {
			sB[i] = (sY[i+1]-sY[i])/sD[i] - sD[i]*(sC[i+1]+2.0*sC[i]);
			sD[i] = (sC[i+1]-sC[i])/sD[i];
			sC[i] = 3.0*sC[i];
    	} sC[sN] = 3.0*sC[sN]; sD[sN] = sD[nm1];

    	for (int p=0;p<sN-1;p++) { sX[p] = sX[p+1]; sY[p] = sY[p+1]; sB[p] = sB[p+1]; sC[p] = sC[p+1]; sD[p] = sD[p+1]; }
    	sX[sN] = 0; sY[sN] = 0; sB[sN] = 0; sC[sN] = 0; sD[sN] = 0;

  		for(l = 0; l < sNN; l++)  retArr[l] = sV[l];
    	i = 0;
    	for(l = 0; l < sNN; l++) {
			ul = retArr[l];
			if(ul < sX[i] || (i < n_1 && sX[i+1] < ul)) {
	    		i = 0; j = sN;
	   			do {
					k = (i+j)/2;
					if(ul < sX[k]) j = k; else i = k; 
	    		} while(j > i+1);
			}
		dx = ul - sX[i]; tmp = (ul < sX[0]) ? 0.0 : sD[i];
		retArr[l] = sY[i] + dx*(sB[i] + dx*(sC[i] + dx*tmp));
    	}
    	
    return(retArr);
	}

    // find min adjustment for all outliers
    public static double[] INTER(double[] inArr, double[] Oarr, double[] fitArr) {

       double[] retArr = new double[inArr.length]; double[] sOarr = new double[Oarr.length]; double[] sfitArr = new double[fitArr.length];
       for (int i=0;i<sOarr.length;i++) { sOarr[i] = Oarr[i]; sfitArr[i] = fitArr[i]; }
       Arrays.sort(sOarr); Arrays.sort(sfitArr);

        for (int i=0;i<inArr.length;i++) {
            int pin = Isearch(sOarr,inArr[i]);
            if (pin < sOarr.length && pin < sfitArr.length) {
                double flo = inArr[i] - sOarr[pin] + sfitArr[pin]; double cel = inArr[i] - sOarr[pin+1] + sfitArr[pin+1]; retArr[i] = (flo + cel) / 2;
            } else {
                retArr[i] = inArr[i] - sOarr[pin] + sfitArr[pin];
            }
        }
        inArr = null; Oarr = null; fitArr = null; sfitArr = null; sOarr = null; System.gc();
    return (retArr);
    }

  	public static void main(String[] args) { 		
		String file = ""; String tfile=""; String format="";
		for (int u=0;u<args.length; u++) { 
		int pin = u+1; String s = args[u]; 
		if (s.equals("-f")) file = args[pin];  
		if (s.equals("-t")) tfile = args[pin];  
		if (s.equals("-fo")) format = args[pin];
		if (s.equals("-th")) thes = Double.parseDouble(args[pin]);
		if (s.equals("-fa")) fact = Double.parseDouble(args[pin]);
		if (s.equals("-kn")) knots = Integer.parseInt(args[pin]);
		if (s.equals("-it")) iter = Integer.parseInt(args[pin]);
		}				
		
		if(format.equals("bed")) {
			double[][] my = null;  my = BED(file); int size = my.length; double[] rx = new double[size]; double[] gy = new double[size]; int nu = size;
    		for (int k=0; k<size;k++) { rx[k] = my[k][0]; gy[k] = my[k][1]; } my = null; System.gc();   	
    		double[] red = eINTER(rx, gy); rx = null; gy = null; System.gc();
        	my = BED(file); rx = new double[size]; gy = new double[size]; for (int k=0; k<size;k++) { rx[k] = my[k][0]; gy[k] = my[k][1]; } my = null; System.gc();
     		double[] green = eINTER(gy, rx); //OutputBED(file, tfile, rx, gy);
     		OutputBED(file, tfile, red, green);
     	}
     	if(format.equals("fe")) {
			double[][] my = null;  my = FE(file); int size = my.length; double[] rx = new double[size]; double[] gy = new double[size]; int nu = size;
    		for (int k=0; k<size;k++) { rx[k] = my[k][0]; gy[k] = my[k][1]; } my = null; System.gc();   	
    		double[] red = eINTER(rx, gy); rx = null; gy = null; System.gc();
        	my = FE(file); rx = new double[size]; gy = new double[size]; for (int k=0; k<size;k++) { rx[k] = my[k][0]; gy[k] = my[k][1]; } my = null; System.gc();
     		double[] green = eINTER(gy, rx); //Output(file, tfile, rx, gy); 
     		Output(file, tfile, red, green);
     	}
  
	}

}
