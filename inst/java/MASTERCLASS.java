/*
 * Author:	Tomas William Fitzgerald
 * Email:	tf2@sanger.ac.uk
 *	All rights are reserved.
 */

import java.io.*;
import java.util.*;
import java.lang.Math.*;
import java.util.ArrayList;
import Sort.*;


class featureObjects{

	public HashMap<String, featureObject> featureHash;
	
	public featureObjects() { featureHash = new HashMap<String, featureObject>(); }
	
	public int size() { return(this.featureHash.size()); }
	
	public void addItem(String k, featureObject FOb) { this.featureHash.put(k, FOb); }

	public featureObject getItem(String k) { return(this.featureHash.get(k)); }
	
	public void removeItem(String k) { this.featureHash.remove(k); }
	
	public Iterator getAllKeys() { return( (Iterator) featureHash.keySet().iterator() ); }

}

class featureObject{

	public String key;
	public double mean, median = 0;
	public int chr, start, stop, len, alnum = 0;
	public Vector indexes, breakPoints, TbreakPoints, algorithms;
	
	public featureObject(String key, Vector breakPoints, Vector algorithms) {
		this.key = key; this.addBreaks(breakPoints); this.addAlgo(algorithms); 
		this.makeUniqueBreaks(); this.makeAlgoNumber();
	}

	public void addBreaks(Vector breakPoints) {
		this.breakPoints=new Vector(); for(int i=0;i<breakPoints.size();i++) { this.breakPoints.add(breakPoints.get(i)); }
	}
	
	public void addAlgo(Vector algorithms) {
		this.algorithms=new Vector(); for(int i=0;i<algorithms.size();i++) { this.algorithms.add(algorithms.get(i)); }
	}
	
	public void makeAlgoNumber() {
		Object[] bps = (Object[]) this.algorithms.toArray(); 
		ArrayList uni = new ArrayList(); java.util.List lis = Arrays.asList(bps);
		for (Object x : lis) { if (!uni.contains(x)) { uni.add(x); } } alnum = uni.size();	
	}
	
	public void makeUniqueBreaks() {
		this.TbreakPoints = new Vector();
		Object[] bps = (Object[]) this.breakPoints.toArray(); 
		ArrayList uni = new ArrayList();
 		java.util.List lis = Arrays.asList(bps);
		for (Object x : lis) {
    		if (!uni.contains(x)) {
        		uni.add(x);  		
    		}
   		}
   		for(int i=0;i<uni.size();i++) {
   			TbreakPoints.add(uni.get(i));
   		} 		
	}
	
	public void makeIndexes(DATAObject data) {
		this.indexes=new Vector();
		Map<String, ArrayList[]> DDATA = data.DATA; 
		for(int i=0;i<this.TbreakPoints.size();i++) {
			//System.out.println(this.key + "\t" + Integer.parseInt(""+TbreakPoints.get(i)));
			int ind = data.findINDEX(this.key, Integer.parseInt(""+TbreakPoints.get(i)));	
   			if(ind < DDATA.get(this.key)[1].size()-1 & i==this.TbreakPoints.size()-1) {
   				while(Integer.parseInt(""+DDATA.get(this.key)[1].get(ind)) == Integer.parseInt(""+DDATA.get(this.key)[1].get(ind+1))) {
   					ind++;
   					//System.out.println(DDATA.get(this.key)[1].get(ind) + "\t" + ind + "\t" + DDATA.get(this.key)[1].size());
   					if(ind >= DDATA.get(this.key)[1].size()-1) { break; } 
   				}
   			}
   			this.indexes.add(ind);
   		} 	
	}
	
	public Vector cutOut(Vector collect, double absCut) {
		
		Vector newDat = new Vector();
		
		if(collect.size() >0) {
			for(int k=0;k<collect.size();k++) {			
				String r1 = ""+collect.get(k); String[] row1 =r1.split("\t");
				double rat = Double.parseDouble(""+row1[3]);
				if(Math.abs(rat) > absCut) {
					newDat.add(collect.get(k));		
				}
			}
		}
		
		return(newDat);
	}

	
	public Vector SplitBreaks(DATAObject data, double absCut) {

		Vector lineOut = new Vector();
		Map<String, ArrayList[]> DATA = data.DATA; 
	
			int st = Integer.parseInt(""+indexes.get(0));
			int so = Integer.parseInt(""+indexes.get(indexes.size()-1));
			//System.out.println(st + "\t" + so);
			//System.out.println(this.TbreakPoints.get(0) + "\t" + this.TbreakPoints.get(TbreakPoints.size()-1));
			if(so>st) {
			Object[] ddata = new Object[(so-st)+1];

			int pin=0;
			for(int j=st;j<=so;j++) {
				ddata[pin] = Double.parseDouble(""+DATA.get(this.key)[2].get(j));
				pin++;
			}
			double M = exMath.mean(ddata);
			double MM = exMath.median(ddata);

		for(int i=0;i<indexes.size()-1;i++) {
			st = Integer.parseInt(""+indexes.get(i));
			so = Integer.parseInt(""+indexes.get(i+1));
		//	System.out.println(st + "\t" + so);
			if(so>st) {

			Object[] rdata = new Object[(so-st)+1];
			pin=0;
			for(int j=st;j<=so;j++) {
				rdata[pin] = Double.parseDouble(""+DATA.get(this.key)[2].get(j));
				pin++;
			}
			
			double m = exMath.mean(rdata);
			double mm = exMath.median(rdata);
//System.out.println(m + "\t" + mm);
				if (Math.abs(m) > absCut || Math.abs(mm) > absCut) {
				//if (Math.abs(m) > absCut) {
					int s1 = st;
					int p = st;
					if(s1 >0) {
						do {
							p = s1;
							if(s1 == 0) { break; }
							s1 = xxReverse(s1, Double.parseDouble(""+DATA.get(this.key)[2].get(s1-1)), MM, absCut);
						} while(s1!=p);
						
					st = s1;
					p = st;
						do {
							p = s1;
							if(s1 == DATA.get(this.key)[0].size()-1) { break; }
							s1 = yyForward(s1, Double.parseDouble(""+DATA.get(this.key)[2].get(s1+1)), MM, absCut);
						} while(s1!=p);
					
					}
					
					int s2 = so;
					p = so;
					if(s2 < DATA.get(this.key)[0].size()) {
						do {
							p = s2;
							if(s2 == DATA.get(this.key)[0].size()-1) { break; }
							s2 = xxForward(s2, Double.parseDouble(""+DATA.get(this.key)[2].get(s2+1)), MM, absCut);
						}  while(s2!=p);
						
					so = s2;
					p = so;
						do {
							p = s2;
							if(s2 == 0) { break; }
							s2 = yyReverse(s2, Double.parseDouble(""+DATA.get(this.key)[2].get(s2-1)), MM, absCut);					
						}  while(s2!=p);
					}
					
					if (s2>s1) {
						Object[] ndata = new Object[(s2-s1)+1];
						pin=0;
						for(int j=s1;j<=s2;j++) {
							ndata[pin] = Double.parseDouble(""+DATA.get(this.key)[2].get(j));
							pin++;
						}
						double new_m = exMath.mean(ndata);
						double new_mm = exMath.median(ndata);
						//System.out.println(this.key + "\t" +  DATA.get(this.key)[0].get(s1) + "\t" + DATA.get(this.key)[1].get(s2)  + "\t" + new_mm + "\t" + ndata.length + "\t" + s1 + "\t" + s2);
					
						lineOut.add(this.key + "\t" +  DATA.get(this.key)[0].get(s1) + "\t" + DATA.get(this.key)[1].get(s2)  + "\t" + new_mm + "\t" + ndata.length + "\t" + s1 + "\t" + s2);
					}

				}

			}
			
		}
		}

		return(lineOut);
	}
	
	public int xxForward(int ind, double rat, double med, double absCut) {
		if(rat > absCut && rat > med / 2) {
			return(ind+1);
		} else if (rat < -absCut && rat < med / 2) {
			return(ind+1);
		} else {
			return(ind);
		}
	}
	
	public int xxReverse(int ind, double rat, double med, double absCut) {
		if(rat > absCut && rat > med / 2) {
			return(ind-1);
		} else if (rat < -absCut && rat < med / 2) {
			return(ind-1);
		} else {
			return(ind);
		}
	}
	
	public int yyForward(int ind, double rat, double med, double absCut) {
		if(rat > absCut && rat > med / 2) {
			return(ind);
		} else if (rat < -absCut && rat < med / 2) {
			return(ind);
		} else {
			return(ind+1);
		}
	}
	
	public int yyReverse(int ind, double rat, double med, double absCut) {
		if(rat > absCut && rat > med / 2) {
			return(ind);
		} else if (rat < -absCut && rat < med / 2) {
			return(ind);
		} else {
			return(ind-1);
		}
	}
	
}

class REPORTObject {

	public static int fieldNum = 5;
	public static ArrayList[] data = null;
	public static ArrayList features = null;

	public REPORTObject(String reportfile) {
		this.getReports(reportfile);
		this.PotentialVariation();
	}
	
	public REPORTObject(Vector result) {
		data = new ArrayList[fieldNum];
		for (int i=0;i<data.length;i++) {
			data[i] = new ArrayList();
		}
		
		for (int i=0;i<result.size();i++) {
        	Vector v = (Vector) result.get(i);

				for(int j=0;j<v.size();j++) {
				    String str = ""+v.get(j);
        			String[] s = str.split("\t");
					for(int k=0;k<data.length;k++) {
						data[k].add(s[k]);
					}
				}		
      	}
      		
      	this.SortAndReplace();
      	this.PotentialVariation();
		
	}

	public static Vector SortTheVector(Vector v1) {
	
		String r1 = ""+v1.get(0); String[] row1 = r1.split("\t");
		String[][] sortData = new String[v1.size()][row1.length];
			
			for(int i=0;i<v1.size();i++) {
				r1 = ""+v1.get(i); row1 = r1.split("\t");
				for(int j=0;j<sortData[0].length;j++) {
					sortData[i][j] = ""+row1[j];
				}
			}
		
		sortData = Sort.SortManager.sort(sortData, new int[]{0, 1}, new int[]{Sort.SortManager.ASC, Sort.SortManager.ASC});	
		
		v1 = new Vector();
		for(int i=0;i<sortData.length;i++) {
			String ro1 = sortData[i][0];
			for(int j=1;j<sortData[0].length;j++) {
				ro1 += "\t" + sortData[i][j];			
			}
			v1.add(ro1);
		}
	return(v1);
	}

	private static void SortAndReplace() {
	
		String[][] sortData = new String[data[0].size()][fieldNum];
	
		for(int i=0;i<data[0].size();i++) {
			for(int j=0;j<fieldNum;j++) {
				sortData[i][j] = ""+data[j].get(i);
			}
		}
		
		sortData = Sort.SortManager.sort(sortData, new int[]{0, 1}, new int[]{Sort.SortManager.ASC, Sort.SortManager.ASC});	
		
		data = new ArrayList[fieldNum];
		for (int i=0;i<fieldNum;i++) {
			data[i] = new ArrayList();
		}
	
		for(int i=0;i<sortData.length;i++) {
			for(int j=0;j<fieldNum;j++) {
				data[j].add(sortData[i][j]);
			}
		}
	
	}
	
	private static ArrayList[] getReports(String filename) {
     
     	data = new ArrayList[fieldNum];
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
        		for(int i=0;i<fieldNum;i++) {
			 		data[i].add(Ldata[i]);
				}			
			} 
				
     		} catch (FileNotFoundException ex) {
         		System.out.println("File not found: " + filename);
     		} catch (IOException ex) {
         		System.out.println(ex.getMessage());
     		} finally {
         		try {
             		if (infile != null) infile.close();
         		} catch (IOException ex) {
             		System.out.println(ex.getMessage());
         		}
     		}
		
		SortAndReplace();
		
      return data;
    }

	// This creates initial feature objects - regions of potential variation. 
	private static ArrayList PotentialVariation() {
	
		features = new ArrayList();
	
		while(data[0].size()>0) {
		
			Object[] ob =  mapBreakPoints(0);	
			features.add(ob);
		}
		return(features);
	}


	// This is an efficent full overlap method - taking items of the pile - data structure must be sorted first!!!
	// - creates non-overlapping feature objects (min start max end and all breakpoint positions)
	private static Vector[] mapBreakPoints(int pin) {
		
		Vector[] values =new Vector[fieldNum];
		for(int i=0;i<fieldNum;i++) {
			values[i] = new Vector();
			values[i].add(data[i].get(0));
			data[i].remove(0);
		}

		for(int i=0;i<values[0].size();i++) {		
		
			boolean go = true;
		
			if(data[0].size()>0) {	
		
				while(go) {

					int chr1 = Integer.parseInt("" + data[0].get(0));
					int chr2 = Integer.parseInt("" + values[0].get(i));
			
					int start1 = Integer.parseInt("" + data[1].get(0));
					int start2 = Integer.parseInt("" + values[1].get(i));
			
					int stop1 = Integer.parseInt("" + data[2].get(0));
					int stop2 = Integer.parseInt("" + values[2].get(i));

					if ( chr1 == chr2 & start1 >= start2 & start1 <= stop2
          				| chr1 == chr2 & stop1 <= stop2 & start1 >= start2
          				| chr1 == chr2 & start1 >= start2 & stop1 <= stop2
          				| chr1 == chr2 & start1 <= start2 & stop1 >= stop2
          				| chr1 == chr2 & start1 <= start2 & stop1 >= start2 ) {
				 				 
				 		for(int j=0;j<data.length;j++) {
				 			values[j].add(data[j].get(0));
				 			data[j].remove(0);
				 		}
				 
				 		if(data[0].size()==0) { go = false; }
				 
					} else {
						go = false;
					}
				}
			}	

		}
		
		return(new Vector[]{values[0], values[1], values[2], values[4]} );
	}		

}


class DATAObject {

	private static int FIELDNUM = 4;
	private static ArrayList[] CHRS = null;
	public static Map<String, ArrayList[]> DATA = null; 

	public DATAObject(String filename) {
		this.getData(filename);
		this.intiDATA();
	}

	public void printDATA() {
		Object[] keys = DATA.keySet().toArray();
   		for(int i=0;i<keys.length;i++) {
   			for(int j =0;j<DATA.get(keys[i])[0].size();j++) {
   				for(int k=0;k<DATA.get(keys[i]).length;k++) {
   					System.out.print(DATA.get(keys[i])[k].get(j) + "\t"); 
   				}
   				System.out.println();
   			}
   		}
	}
	
	public int findINDEX(String key, int start) {
		
		int index = 0;
		boolean found = false;
		int i=0;

		while(!found) {
			int cstart = Integer.parseInt("" + DATA.get(key)[0].get(i) );
			if(cstart == start) { found = true; index = i; }
			int cstop = Integer.parseInt("" + DATA.get(key)[1].get(i) );
			if(cstop == start) { found = true; index = i; }			
			i++;
			if( i==DATA.get(key)[0].size() ) {
				int pin=i-1; int check1 = start-Integer.parseInt("" + DATA.get(key)[1].get(pin)); pin--;
				while(check1<0) { check1 = start-Integer.parseInt("" + DATA.get(key)[1].get(pin)); pin--; }
				index=pin; break;
			}
		}
		
	return(index);
	}
    
    private void intiDATA() {
         
		DATA = new HashMap<String, ArrayList[]>();
		
    	Object[] chrs = new Object[CHRS[0].size()];
 		for(int i=0;i<chrs.length;i++) {
 			chrs[i] = CHRS[0].get(i); 
 		}
         
        ArrayList uni = new ArrayList();
        java.util.List lis = Arrays.asList(chrs);
		for (Object x : lis) {
    		if (!uni.contains(x)) {
        		uni.add(x);
    		}
   		}
        
        for(int i=0;i<uni.size();i++) {
        	ArrayList[] thisData = new ArrayList[FIELDNUM-1];
        	for(int x=0;x<FIELDNUM-1;x++) {
        		thisData[x] = new ArrayList();
        	}
        	for(int j=0;j<CHRS[0].size();j++) {
        		if(CHRS[0].get(j).equals(uni.get(i))) {
        			int pin = 0;
        			for(int k=1;k<CHRS.length;k++) {
        				thisData[pin].add(CHRS[k].get(j));
        				pin++;
        			}
        		}
        	}
        DATA.put(""+uni.get(i), thisData);
        } 
    CHRS=null; System.gc();
    }
    
	private void getData(String filename) {

     	filename.trim();
     	String inLine;
     	BufferedReader infile = null;
     
    	CHRS = new ArrayList[FIELDNUM];
     	for(int i=0;i<FIELDNUM;i++) {
     		CHRS[i] = new ArrayList();
     	}
     
     try {
     
         infile = new BufferedReader(new FileReader (filename));
         while ((inLine=infile.readLine()) != null) {   	
        	String[] Ldata = inLine.split("\t");
			for(int j=0;j<FIELDNUM;j++) {
				CHRS[j].add(Ldata[j]);;
			}
         }
         
        infile.close();
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
  }
}


class MASTERCLASS {

	private static DATAObject DATA;
	private static REPORTObject REPORT;
	private static String outfile;
	private static double mParam=0.5;
	private static double absCut=0.3;
	private static double s_param = 2.2;
	private static double r_param = 0.5;
	
	public MASTERCLASS(String datafile, String reportfile, String outfile) {
		this.DATA=new DATAObject(datafile);
		this.REPORT=new REPORTObject(reportfile);
		this.outfile=outfile;
	}
	
	public MASTERCLASS(String datafile, String reportfile, String outfile, double abscut) {
		this.absCut=absCut;
		this.DATA=new DATAObject(datafile);
		this.REPORT=new REPORTObject(reportfile);
		this.outfile=outfile;
	}
	
	public MASTERCLASS(String datafile, String reportfile, String outfile, double mParam, double abscut, double s_param, double r_param) {
		this.DATA=new DATAObject(datafile);
		this.REPORT=new REPORTObject(reportfile);
		this.outfile=outfile;
		this.mParam=mParam;
		this.absCut=absCut;
		this.s_param=s_param;
		this.r_param=r_param;
	}
	
	private void Out(Vector result) {
        
        PrintWriter pw = null;
    	File output = new File(outfile);    

    	try {
        	pw = new PrintWriter(new FileOutputStream(output), true);
        	
        	for (int i=0;i<result.size();i++) {
				pw.print(result.get(i) + "\n");	
      		}
    	}
    	catch (IOException ex) {
        	System.out.println(ex.getMessage());
    	}
    	finally {
        	if (pw != null) pw.close();
    	}
    }
	
	private Vector mSegments() {
	
		System.out.println(absCut);

		Vector finV = new Vector();
		ArrayList features = REPORT.features;
		
			for(int i=0;i<features.size();i++) {	
		
				Vector[] v=  (Vector[]) features.get(i);		
				String key = "" + v[0].get(0);
			//System.out.println(key);
					Vector breakps = new Vector();
					for( int j = 0;j<v[0].size();j++) {
						breakps.add(v[1].get(j));
						breakps.add(v[2].get(j));
					}
					
					int[] bps = new int[breakps.size()];
					for(int f=0;f<bps.length;f++) {
						bps[f] = Integer.parseInt("" + breakps.get(f));
					}					
					Arrays.sort(bps);
					Vector breaks = new Vector();
					for(int j=0;j<bps.length;j++) {
						breaks.add(bps[j]);
					}
					Vector algo = new Vector();
					for(int j=0;j<v[3].size();j++) {
						algo.add(v[3].get(j));
					}

					featureObject fOb = new featureObject(key, breaks, algo);
					fOb.makeIndexes(this.DATA);
					finV.add(fOb.cutOut(fOb.SplitBreaks(this.DATA, this.absCut), this.absCut));
			}
			
		return(finV);	
	}
    
    private Vector mCondense(REPORTObject report) {
	
		Vector finV = new Vector();
		ArrayList features = report.features;
		Map<String, ArrayList[]> DDATA = this.DATA.DATA; 
		
		for(int i=0;i<features.size();i++) {	
		
				Vector[] v=  (Vector[]) features.get(i);			
				String key = "" + v[0].get(0);
			
					Vector breakps = new Vector();
					for( int j = 0;j<v[0].size();j++) {
						breakps.add(v[1].get(j));
						breakps.add(v[2].get(j));			
					}
					
					int[] bps = new int[breakps.size()];
					for(int f=0;f<bps.length;f++) {
						bps[f] = Integer.parseInt("" + breakps.get(f));
					}

					Arrays.sort(bps);
					Vector breaks = new Vector();
					for(int j=0;j<bps.length;j++) {
						breaks.add(bps[j]);
					}
					Vector algo = new Vector();
					for(int j=0;j<v[3].size();j++) {
						algo.add(v[3].get(j));
					}

					featureObject fOb = new featureObject(key, breaks, algo);
					fOb.makeIndexes(this.DATA);
	
					int st = Integer.parseInt(""+fOb.indexes.get(0));
					int so = Integer.parseInt(""+fOb.indexes.get(fOb.indexes.size()-1));
					Object[] ddata = new Object[(so-st)+1];

					int pin=0;
					for(int j=st;j<=so;j++) {
						ddata[pin] = Double.parseDouble(""+DDATA.get(key)[2].get(j));
						pin++;
					}
					double MM=0;
					if(pin>=1) {
						MM = exMath.median(ddata);
					}
					finV.add(key + "\t" + DDATA.get(key)[0].get(st) + "\t" +  DDATA.get(key)[1].get(so) + "\t" +  MM + "\t" +  ddata.length + "\t" + st + "\t" + so);
				}

		return(finV);
	}
	
	private boolean mermer(String[] row1, String[] row2, int s1) {
			
			boolean merge = false;
			double m1 = Double.parseDouble(""+row1[3]);
			double m2 = Double.parseDouble(""+row2[3]);
			double size = Math.sqrt(Integer.parseInt(""+row1[4]))*s_param;
			
			if(m1>0) {
				if(m2>(m1*r_param) && size > s1) {			
					merge = true;
				}
			} else if (m1<0) {
				if(m2<(m1*r_param) && size > s1) {
					merge = true;
				} 
			}
		return(merge);		
	}
	
	private Vector mer(Vector v1) {
	
		Map<String, ArrayList[]> DDATA = DATA.DATA; 
		Vector v = new Vector();
		Vector v2 = new Vector();
		
		for(int i=0;i<v1.size()-1;i++) {
		
			String r1 = ""+v1.get(i); String[] row1 = r1.split("\t");
			String r2 = ""+v1.get(i+1); String[] row2 = r2.split("\t");	
			int s1 = (Integer.parseInt(""+row2[5]))-(Integer.parseInt(""+row1[6]));
			
			if(mermer(row1, row2, s1)) {
					int st = Integer.parseInt(""+row1[5]);
					int so = Integer.parseInt(""+row2[6]);
					//System.out.println(row1[0] + "\t" + st + "\t" + so);
				if(so>st) {
					Object[] rata = new Object[(so-st)+1];			
					int pin=0;
					for(int j=st;j<=so;j++) {
						rata[pin] = Double.parseDouble(""+DDATA.get(row1[0])[2].get(j));
						pin++;
					}
					double MM = exMath.median(rata);
			
				v2.add(row1[0] + "\t" + row1[1] + "\t" + row2[2] + "\t" + MM + "\t" + rata.length + "\t" + row1[5] + "\t" + row2[6]);
				}
			} else {
				v2.add(r1);
				v2.add(r2);
			}
		}
		if(v1.size()==1) { v2.add(v1.get(0)); }
		v2 = REPORTObject.SortTheVector(v2);
		v.add(v2); 
		v2 = mCondense(new REPORTObject(v));
		Vector v3 = new Vector();
		
		for(int i=v2.size()-1;i>0;i--) {
		
			String r1 = ""+v2.get(i); String[] row1 = r1.split("\t");
			String r2 = ""+v2.get(i-1); String[] row2 = r2.split("\t");
			int s1 = (Integer.parseInt(""+row1[5]))-(Integer.parseInt(""+row2[6]));

			if(mermer(row1, row2, s1)) {
					int st = Integer.parseInt(""+row2[5]);
					int so = Integer.parseInt(""+row1[6]);
				if(so>st) {
					Object[] rata = new Object[(so-st)+1];
					int pin=0;
					for(int j=st;j<=so;j++) {
						rata[pin] = Double.parseDouble(""+DDATA.get(row1[0])[2].get(j));
						pin++;
					}
					double MM = exMath.median(rata);
			
				v3.add(row1[0] + "\t" + row2[1] + "\t" + row1[2] + "\t" + MM + "\t" + rata.length + "\t" + row2[5] + "\t" + row1[6]);
				}
			} else {
				v3.add(r1);
				v3.add(r2);
			}
		}
		if(v2.size()==1) { v3.add(v2.get(0)); }
		v3 = REPORTObject.SortTheVector(v3);
		v = new Vector();
		v.add(v3);
		
		return(mCondense(new REPORTObject(v)));
	}

	private void Solidate() {
	
		Vector rr=new Vector();		
		Vector result1 = mSegments();

		Vector result=new Vector();
		for (int i=0;i<result1.size();i++) {
        	Vector v = (Vector) result1.get(i);
				for(int j=0;j<v.size();j++) {
				    String str = ""+v.get(j);
        			String[] s = str.split("\t");
					result.add(str);
				}		
      	}
			if(result.size()>0) {
				result = REPORTObject.SortTheVector(result);
				Vector v1 = new Vector(); v1.add(result);
				REPORTObject report = new REPORTObject(v1);
				result1 = mCondense(report);
				do {		
					result = result1;
					result1 = mer(result1);
					//System.out.println(result.size() + "\t" + result1.size());
				} while(result.size() > result1.size());
			Out(result1);
			} else {
				Out(new Vector());
			}
	}


	public static void main(String[] args) {
	
		String datafile = "";
		String reportfile = "";	
		String outfile = "";
		double abscut=0;
		
		for (int u=0;u<args.length; u++) {		
			int pin = u+1;
            String s = args[u]; 			
			if (s.equals("-f")) { 
				datafile = args[pin];
			}
			if (s.equals("-r")) { 
				reportfile = args[pin];
			} 
			if (s.equals("-o")) { 
				outfile = args[pin];
			} 
			if (s.equals("-abs")) { 
				absCut = Double.parseDouble(""+args[pin]);
			} 
			 
		}	
		MASTERCLASS m = new MASTERCLASS(datafile, reportfile, outfile, abscut);
		m.Solidate();
		//Vector result = Solidate(datafile, reportfile);
		//Out(outfile, result);
	}


}







