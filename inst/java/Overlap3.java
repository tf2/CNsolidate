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

public class Overlap3 {
		
	public static String Gold_standard = "42Mcalls_all_hg19.txt";
	
	/* Default constructor */
	public Overlap3() {

	}

	    public static int OverOverLap(int chr1, int chr2, int start1, int start2, int stop1, int stop2) {
        int check;
        
        if (chr1 == chr2 & start1 >= start2 & start1 <= stop2
          | chr1 == chr2 & stop1 <= stop2 & start1 >= start2
          | chr1 == chr2 & start1 >= start2 & stop1 <= stop2
          | chr1 == chr2 & start1 <= start2 & stop1 >= stop2
          | chr1 == chr2 & start1 <= start2 & stop1 >= start2) {
            
            check = 1;
        } 
            else {
            
            check =0;
        }

        return check;
    }


 public static String[][] ProcessAll(int[] chr1, int[] chr2, int[] start1, int[] start2, int[] stop1, int[] stop2, String[] features) {
    
        String[][] over = new String[chr1.length][3];   
        String overlapping_feature = "";    
            for (int i=0;i<chr1.length ; i++) {
                int pin = 0;
				double dist = 0;
				double rat = 0;
				//overlapping_feature = features[i];
				
             for (int k =0 ;k<chr2.length ; k++) {
                 int pin2 = OverOverLap(chr1[i], chr2[k], start1[i], start2[k], stop1[i], stop2[k]);
                 	                 	
                 	              
                 	if (pin2 == 1) { 
                 		pin += pin2;
                 	
                 		//System.out.println(chr1[i] + "\t" + chr2[k] + "\t" + start1[i] + "\t" + start2[k] + "\t" + stop1[i] + "\t" + stop2[k] + "\t" + pin);
                 		if (start1[i] <= start2[k] && stop1[i] >= stop2[k]) {
                 		
                 			double con_dist = stop2[k] - start2[k]; 
                 			double test_dist = stop1[i] - start1[i];
                 			double t = (con_dist/test_dist)*100;
                 			dist += t; 
                 			
                 		} else if (start1[i] >= start2[k] && stop1[i] <= stop2[k]) {
                 			dist += 100;
                 		} else if (start1[i] <= start2[k] && stop1[i] >= start2[k] && stop1[i] <= stop2[k]) {
                 		      		
                 			double con_dist = stop1[i] - start2[k]; 
                 			double test_dist = stop1[i] - start1[i]; 
                 		    double t = (con_dist/test_dist)*100;
                 			dist += t;
                 			
                 		} else if (stop1[i] >= stop2[k] && start1[i] <= stop2[k] && stop1[i] >= start2[k]) {        		
                 			double con_dist = stop2[k] - start1[i]; 
                 			double test_dist = stop1[i] - start1[i]; 
                 		    double t = (con_dist/test_dist)*100;
                 			dist += t;
                 		}
                 	} 
             }
                if (pin > 0) {                 	
                	if (dist > 100) { dist = 100; }
                	over[i][0] = ""+pin; 
                	over[i][1] = ""+dist;
                	over[i][2] = "";//+overlapping_feature;
                }
                else if (pin == 0) {              
                	over[i][0] = ""+0; 
                	over[i][1] = ""+0;
                	over[i][2] = ""+0;
                }
            }
      
        return over;
    }

     
    public static int[] MakeFDR(String filename, String outname, String type) {
        
        	int[] inter = new int[2]; inter[0]=0; inter[1]=0;
			String FDR="";
        	
        	int col = testNumberCol(filename,1);
        	ArrayList[] Data = getData(filename, col, 1);
        	
        	//System.out.println("hi");
        	ArrayList[] Gold = getData(Gold_standard, 3, 1);
  		
  			int[] chr1 = new int[Gold[0].size()];
			int[] start1 = new int[Gold[0].size()];
			int[] stop1 = new int[Gold[0].size()];
			String[] men1 = new String[Gold[0].size()];

			for (int i=0;i<chr1.length;i++) {
				//System.out.println(Gold[0].get(i) + "\t" + Gold[1].get(i) + "\t" + Gold[2].get(i));
				chr1[i] = Integer.parseInt("" + Gold[0].get(i));
				start1[i] = Integer.parseInt("" + Gold[1].get(i));
				stop1[i] = Integer.parseInt("" + Gold[2].get(i));
				men1[i] = "" + Gold[2].get(i);
			}
				
				int[] chr2 = new int[Data[0].size()];
  				int[] start2 = new int[Data[0].size()];
  				int[] stop2 = new int[Data[0].size()];
  				//String[] men2 = new String[Data[0].size()];
  			
  				for (int j=0;j<chr2.length;j++) {
        			chr2[j] = Integer.parseInt("" + Data[0].get(j));
					start2[j] = Integer.parseInt("" + Data[1].get(j));
					stop2[j] = Integer.parseInt("" + Data[2].get(j));
					//men2[j] = Double.parseDouble("" + 0);
  				}
  				
  				//System.out.println("check1");
  				
				String[][] over = ProcessAll(chr2, chr1, start2, start1, stop2, stop1,men1);
				//FDR += Clean.MakeTable("False_"+algo, over, "");
				//System.out.println(filename);
				String[] values = new String[over.length];
				ArrayList outt = new ArrayList();
				for(int i=0;i<over.length;i++) {
					
					String row = "";
					for(int k=0; k<col-1;k++) {
						row = row + Data[k].get(i) + "\t";
					}
					row = row + Data[col-1].get(i) + "\t" +  over[i][1]  + "\n";
					outt.add(row);
					//System.out.print(row);
					
					values[i]=row;
					if(Integer.parseInt(""+over[i][0])==0) { inter[0]++; } 
					else { inter[1]++; }
				}
				Handles hand = new Handles(filename, outname);
				hand.ArrayListOut(outt);
				
		return(inter);		
    }
    
    
    
    public static void Output( String filename, String[] values) {
    PrintWriter pw = null; File output = new File(filename); 
    if(output.exists()) System.out.println("File already exists - overwriting!!!");
        try {
        	pw = new PrintWriter(new FileOutputStream(output), true);
 			for(int i=0;i<values.length;i++) {
 				pw.print(values[i]); 
 			}
 		}
        catch (IOException ex) { System.out.println(ex.getMessage());}
    	finally { if (pw != null) pw.close(); }
     }
    
    
    public static int testNumberCol(String filename, int flag) {
     filename.trim();
     String inLine;
     BufferedReader infile = null;
     int count = 0; int rcount=0;
     
     try {
         infile = new BufferedReader(new FileReader (filename));
   
         while ((inLine=infile.readLine()) != null & rcount==0) {
        	String[] Ldata = null;
        	if(flag==0) {
        		Ldata = inLine.split(" ");
        	} else if(flag==1) {
        		Ldata = inLine.split("\t");
        	}
        	for (int i=0;i<Ldata.length;i++) {
        		count++;
			 }
			 rcount++;
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

      return count;

    }
    
    public static ArrayList[] getData(String filename, int ncol, int flag) {
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
        	String[] Ldata = null;
        	if(flag==0) {
        		Ldata = inLine.split(" ");
        	} else if(flag==1) {
        		Ldata = inLine.split("\t");
        	}
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
    
    public static void runDirFDR(String dir, String outname) {
    	int[] p = MakeFDR(dir, outname, "false");
    }
    

	public static void main(String[] args) {
	
		String file= "";
		String out= "";
		String al= "";
		
		for (int u=0;u<args.length; u++) {
			
			int pin = u+1;
            String s = args[u]; 
			
			if (s.equals("-f")) { 
				file = args[pin];
			} 
			if (s.equals("-o")) { 
				out = args[pin];
			} 
			if(s.equals("-g")) { 
				Gold_standard = args[pin];
			}
		}	
		runDirFDR(file, out);

	}


}

