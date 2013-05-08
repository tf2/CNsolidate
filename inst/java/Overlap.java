/**
 *
 *  This program is distributed in the hope that it might be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  http://www.r-project.org/Licenses/
 *
*/

import java.io.*;
import java.util.*;
import java.lang.Math.*;
import java.util.ArrayList;

public class Overlap {
		
	public static String Gold_standard = "../data/42Mcalls_all_feb09.txt";
	
	/* Default constructor */
	public Overlap() {

	}

    public static int OverOverLap(int chr1, int chr2, int start1, int start2, int stop1, int stop2) {
        int check;      
        if (chr1 == chr2 & start1 >= start2 & start1 <= stop2
          | chr1 == chr2 & stop1 <= stop2 & start1 >= start2
          | chr1 == chr2 & start1 >= start2 & stop1 <= stop2
          | chr1 == chr2 & start1 <= start2 & stop1 >= stop2
          | chr1 == chr2 & start1 <= start2 & stop1 >= start2) {     
        check = 1;
        } else { check =0; }
        return check;
    }

    public static double[][] ProcessAll(int[] chr1, int[] chr2, int[] start1, int[] start2, int[] stop1, int[] stop2) {    
        double[][] over = new double[chr1.length][2];            
            for (int i=0;i<chr1.length ; i++) {
                int pin = 0; double dist = 0;			
             for (int k =0 ;k<chr2.length ; k++) {
                 int pin2 = OverOverLap(chr1[i], chr2[k], start1[i], start2[k], stop1[i], stop2[k]);            	
                 	if (pin2 == 1) { 
                 		pin += pin2;
                 		if (start1[i] <= start2[k] && stop1[i] >= stop2[k]) {
                 			double con_dist = stop2[k] - start2[k]; double test_dist = stop1[i] - start1[i];
                 			double t = (con_dist/test_dist)*100; dist += t;             			
                 		} else if (start1[i] >= start2[k] && stop1[i] <= stop2[k]) { dist += 100;
                 		} else if (start1[i] <= start2[k] && stop1[i] >= start2[k] && stop1[i] <= stop2[k]) {        		
                 			double con_dist = stop1[i] - start2[k]; double test_dist = stop1[i] - start1[i]; 
                 		    double t = (con_dist/test_dist)*100; dist += t;          			
                 		} else if (stop1[i] >= stop2[k] && start1[i] <= stop2[k] && stop1[i] >= start2[k]) {           		
                 			double con_dist = stop2[k] - start1[i]; double test_dist = stop1[i] - start1[i]; 
                 		    double t = (con_dist/test_dist)*100; dist += t;
                 		}
                 	} 
             }
                if (pin > 0) {  if (dist > 100) { dist = 100; } over[i][0] = pin; over[i][1] = dist; }
                else if (pin == 0) { over[i][0] = 0; over[i][1] = 0; }
            }
        return over;
    }


 public static double[][] ProcessAll2(int[] chr1, int[] chr2, int[] start1, int[] start2, int[] stop1, int[] stop2, double[] rats) {    
        double[][] over = new double[chr1.length][3];            
            for (int i=0;i<chr1.length ; i++) {
                int pin = 0; double dist = 0;double rat = 0;		
             for (int k =0 ;k<chr2.length ; k++) {
                 int pin2 = OverOverLap(chr1[i], chr2[k], start1[i], start2[k], stop1[i], stop2[k]);                 	                               	              
                 	if (pin2 == 1) { 
                 		pin += pin2; rat += rats[k];
                 		if (start1[i] <= start2[k] && stop1[i] >= stop2[k]) {                	
                 			double con_dist = stop2[k] - start2[k]; double test_dist = stop1[i] - start1[i];
                 			double t = (con_dist/test_dist)*100; dist += t;                			
                 		} else if (start1[i] >= start2[k] && stop1[i] <= stop2[k]) {
                 			dist += 100;
                 		} else if (start1[i] <= start2[k] && stop1[i] >= start2[k] && stop1[i] <= stop2[k]) {                 		      		
                 			double con_dist = stop1[i] - start2[k]; double test_dist = stop1[i] - start1[i]; 
                 		    double t = (con_dist/test_dist)*100; dist += t;                 			
                 		} else if (stop1[i] >= stop2[k] && start1[i] <= stop2[k] && stop1[i] >= start2[k]) {        		
                 			double con_dist = stop2[k] - start1[i]; double test_dist = stop1[i] - start1[i]; 
                 		    double t = (con_dist/test_dist)*100; dist += t;
                 		}
                 	} 
             }
                if (pin > 0) {                 	
                	if (dist > 100) { dist = 100; }
                	over[i][0] = pin; over[i][1] = dist; over[i][2] = rat/pin;
                }
                else if (pin == 0) {              
                	over[i][0] = 0; over[i][1] = 0; over[i][2] = 0;
                }
            }   
        return over;
    }

     
    	public static void MakeFinalReport(String filename, ArrayList uni, ArrayList[] result) {
        
        PrintWriter pw = null;
    	File output = new File(filename);    


    	try {
        	pw = new PrintWriter(new FileOutputStream(output), true);
        	
        	for(int i=0;i<uni.size();i++) {
				if(i<uni.size()-1) {
					pw.print(uni.get(i) + "\t");
				} else {
					pw.print(uni.get(i) + "\n");
				}
			}
			for(int i=0;i<result[0].size();i++) {
				for(int j=0;j<result.length;j++) {
					if(j<result.length-1) {
						pw.print(result[j].get(i) + "\t");
					} else {
						pw.print(result[j].get(i) + "\n");
					}
				}
			}
        	pw.close();
        	
    	}
    	catch (IOException ex) {
        	System.out.println(ex.getMessage());
    	}
    	finally {
        	if (pw != null) pw.close();
    	}
    }
    
    public static void MakeFDR(String filename, ArrayList[] result, int num) {
        

			String FDR="";

        	ArrayList[] Gold = getData(Gold_standard, 3);
  		
  			int[] chr1 = new int[Gold[0].size()];
			int[] start1 = new int[Gold[0].size()];
			int[] stop1 = new int[Gold[0].size()];
			double[] men1 = new double[Gold[0].size()];

			for (int i=0;i<chr1.length;i++) {
				chr1[i] = Integer.parseInt("" + Gold[0].get(i));
				start1[i] = Integer.parseInt("" + Gold[1].get(i));
				stop1[i] = Integer.parseInt("" + Gold[2].get(i));
				men1[i] = Double.parseDouble("" + 0);
			}
        	

			for(int k = 0; k<num;k++) {
			
				int thisN = k+1;
				
				Vector chr = new Vector();
  				Vector sta = new Vector();
  				Vector sto = new Vector();
  				Vector men = new Vector();

				for(int i=0;i<result[0].size();i++) {
					int c =Integer.parseInt("" + result[0].get(i));
					int n = (int)Double.parseDouble(""+result[result.length-1].get(i));
					//System.out.println(n + "\t" + thisN);
					if(n >= thisN &&  c < 23) {
						chr.add(result[0].get(i));
  						sta.add(result[1].get(i));
  						sto.add(result[2].get(i));
  						men.add(result[3].get(i));
					}
				}
				
				int[] chr2 = new int[chr.size()];
  				int[] start2 = new int[chr.size()];
  				int[] stop2 = new int[chr.size()];
  				double[] men2 = new double[chr.size()];
  			
  				for (int j=0;j<chr2.length;j++) {
        			chr2[j] = Integer.parseInt("" + chr.get(j));
					start2[j] = Integer.parseInt("" + sta.get(j));
					stop2[j] = Integer.parseInt("" + sto.get(j));
					men2[j] = Double.parseDouble("" + men.get(j));
  				}
  				
  			//	System.out.println(chr2.length + "\t" + chr1.length);
				double[][] over = ProcessAll2(chr2, chr1, start2, start1, stop2, stop1,men1);
				FDR += Handles.MakeTable(""+thisN, over, "");
				
			}
			
		PrintWriter pw = null;
    	File output = new File(filename);    
    	try {
        	pw = new PrintWriter(new FileOutputStream(output), true);        	
			pw.print(FDR);		      	
        	pw.close();
        	
    	}
    	catch (IOException ex) {
        	System.out.println(ex.getMessage());
    	}
    	finally {
        	if (pw != null) pw.close();
    	}
    }
    
    public static ArrayList[] UpdateReport(String filename, ArrayList[] data, ArrayList[] nums) {
        
        PrintWriter pw = null;
    	File output = new File(filename);    
		ArrayList[] Fdata = new ArrayList[data.length+1];
		for(int i=0;i<Fdata.length;i++) {
			Fdata[i] = new ArrayList();
		}

    	try {
        	pw = new PrintWriter(new FileOutputStream(output), true);
        	

			for(int i=0;i<data[0].size();i++) {
				for(int j=0;j<data.length;j++) {
					if(j<data.length-1) {
						pw.print(data[j].get(i) + "\t");
						//System.out.print(data[j].get(i) + "\t");
						Fdata[j].add(data[j].get(i));
					} else {
						int sum = 0;
						for(int k=0;k<nums.length;k++) {
							sum+=Double.parseDouble(""+nums[k].get(i));
						}
							pw.print(data[j].get(i) + "\t" + sum + "\n");
							//System.out.print(data[j].get(i) + "\t" + sum + "\n");
							Fdata[j].add(data[j].get(i));
							Fdata[j+1].add(sum);		
					}
				}
			}
        	pw.close();
        	
    	}
    	catch (IOException ex) {
        	System.out.println(ex.getMessage());
    	}
    	finally {
        	if (pw != null) pw.close();
    	}
    	return(Fdata);
    }

    public static ArrayList[] getData(String filename, int ncol) {
     filename.trim();
     String inLine;
     BufferedReader infile = null;
     
     	ArrayList[] data = new ArrayList[ncol];
		for (int i=0;i<data.length;i++) {
			data[i] = new ArrayList();
		}
     
     try {
         infile = new BufferedReader(new FileReader (filename));
   
         while ((inLine=infile.readLine()) != null) {
        	
        	String[] Ldata = inLine.split("\t");
        	for (int i=0;i<ncol;i++) {
        		data[i].add(Ldata[i]);
			 }
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
	
	
	public static void Run(String filename1, String filename2, String path) {
		String command = "mkdir " + path + "Cmatrices";
		Clean.ExeShellCommand(command);
		command = "chmod 777 " + path + "Cmatrices";
		Clean.ExeShellCommand(command);
		
		ArrayList[] Adata = getData(path + filename1, 5);
		int ind = filename1.lastIndexOf("_"); 
		String name = filename1.substring(0, ind);
		String repName1 = path + "Cmatrices/" + name + "_NumberMatrix.txt";
		String repName2 = path + "Cmatrices/" + name + "_PercentageMatrix.txt";
		String repName3 = path + "Cmatrices/" + name + "_RatioMatrix.txt";
		String repName4 = path + "Cmatrices/" + name + "_SizeMatrix.txt";
		String repName5 = path + "Cmatrices/" + name + "_FDR.txt";
		
		Vector types = new Vector();
		
		for (int i=0;i<Adata[0].size();i++) {
			types.add(Adata[4].get(i));
  		}
  		
  		ArrayList[] Rdata = getData(path + filename2, 7);
  		
  		int[] chr2 = new int[Rdata[0].size()];
		int[] start2 = new int[Rdata[0].size()];
		int[] stop2 = new int[Rdata[0].size()];
		
		for (int i=0;i<chr2.length;i++) {	
			chr2[i] = Integer.parseInt("" + Rdata[0].get(i));
			start2[i] = Integer.parseInt("" + Rdata[1].get(i));
			stop2[i] = Integer.parseInt("" + Rdata[2].get(i));
		}
  		
  		Object[] uniT = (Object[]) types.toArray(); 
		ArrayList uni = new ArrayList();
 		java.util.List lis = Arrays.asList(uniT);
		for (Object x : lis) {
    		if (!uni.contains(x)) {
        		uni.add(x);        		
    		}
   		}
  		
  		ArrayList[] OneZresult = new ArrayList[uni.size()];
  		ArrayList[] Numresult = new ArrayList[uni.size()];
  		ArrayList[] Ratresult = new ArrayList[uni.size()];
  		ArrayList[] Perresult = new ArrayList[uni.size()];
  		ArrayList[] Sizresult = new ArrayList[uni.size()];
  		for(int i=0;i<OneZresult.length;i++) {
  			OneZresult[i] = new ArrayList();
  			Numresult[i] = new ArrayList();
  			Ratresult[i] = new ArrayList();
  			Perresult[i] = new ArrayList();
  			Sizresult[i] = new ArrayList();
  		}
  		
  		for(int i=0;i<uni.size();i++) {
  			Vector chr = new Vector();
  			Vector sta = new Vector();
  			Vector sto = new Vector();
  			Vector men = new Vector();
  			String check = "" + uni.get(i);
  			for(int j=0;j<types.size();j++) {
  				int c = Integer.parseInt(""+Adata[0].get(j));
  					//if(types.get(j).equals(check) & c < 23) {
  					if(types.get(j).equals(check)) {
  						chr.add(Adata[0].get(j));
  						sta.add(Adata[1].get(j));
  						sto.add(Adata[2].get(j));
  						men.add(Adata[3].get(j));
  					}
  			}
  			
  			int[] chr1 = new int[chr.size()];
  			int[] start1 = new int[chr.size()];
  			int[] stop1 = new int[chr.size()];
  			double[] men1 = new double[chr.size()];
  			
  			for (int j=0;j<chr1.length;j++) {
        		chr1[j] = Integer.parseInt("" + chr.get(j));
				start1[j] = Integer.parseInt("" + sta.get(j));
				stop1[j] = Integer.parseInt("" + sto.get(j));
				men1[j] = Double.parseDouble("" + men.get(j));
  			}
  			
			double[][] over = ProcessAll2(chr2, chr1, start2, start1, stop2, stop1,men1);
			for(int j = 0;j<over.length;j++) {
				//if(chr2[j] <23) {
					double p = over[j][0];
					if(p>1) { p = 1; }
  					OneZresult[i].add(p);
  					Numresult[i].add(over[j][0]);
  					Perresult[i].add(over[j][1]);
  					Ratresult[i].add(over[j][2]);
  					Sizresult[i].add(stop2[j]-start2[j]);
  				//}
  			}
  		}
		
		MakeFinalReport(repName1, uni, Numresult);
		MakeFinalReport(repName2, uni, Perresult);
		MakeFinalReport(repName3, uni, Ratresult);
		MakeFinalReport(repName4, uni, Sizresult);
		MakeFDR(repName5, UpdateReport(path + filename2, Rdata, OneZresult), Numresult.length);
		
	}


	public static void main(String[] args) {
	
		String file1= "";
		String file2= "";
		
		for (int u=0;u<args.length; u++) {
			
			int pin = u+1;
            String s = args[u]; 
			
			if (s.equals("-f1")) { 
				file1 = args[pin];
			} 
			
			if (s.equals("-f2")) { 
				file2 = args[pin];	
			}
			if(s.equals("-g")) { 
				Gold_standard = args[pin];
			}
		}	

		String path="";
		String justFileName1="";
		String justFileName2="";
		int index = file1.lastIndexOf('/');

		if (index >= 0) {
            justFileName1 = file1.substring(index + 1);
            path = file1.substring(0, index) + "/";
        } else {
            index = file1.lastIndexOf('\\');
            if (index >= 0) {
                justFileName1 = file1.substring(index + 1) + "\\";
                path = file1.substring(0, index);
            } else {
                justFileName1 = file1;
                path = file1;
            }
        }
        
        index = file2.lastIndexOf('/');

		if (index >= 0) {
            justFileName2 = file2.substring(index + 1);
        } else {
            index = file2.lastIndexOf('\\');
            if (index >= 0) {
                justFileName2 = file2.substring(index + 1) + "\\";
            } else {
                justFileName2 = file2;
            }
        }

		//System.out.println(path);


		Run(justFileName1, justFileName2, path);
	}


}

