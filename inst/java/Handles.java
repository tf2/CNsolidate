import java.io.*;
import java.util.*;
import java.util.ArrayList;

public class Handles {

	private String inFileName = "";
	private String outFileName = "";
	
	public Handles(String inFileName, String outFileName) {
		this.inFileName=inFileName; this.outFileName=outFileName;		
	}
	
	public ArrayList[] getData(int ncol) {
    	inFileName.trim(); String inLine; BufferedReader infile = null; 
     	ArrayList[] data = new ArrayList[ncol]; for (int i=0;i<data.length;i++) { data[i] = new ArrayList(); }
     
     	try {
         	infile = new BufferedReader(new FileReader (inFileName));  
         	while ((inLine=infile.readLine()) != null) { 	
        		String[] Ldata = inLine.split("\t");
        		for (int i=0;i<ncol;i++) { data[i].add(Ldata[i]); }
         	}
     	} catch (FileNotFoundException ex) { System.out.println("File not found: " + inFileName); }
     	catch (IOException ex) { System.out.println(ex.getMessage()); }
     	finally { try { if (infile != null) infile.close(); }
         	catch (IOException ex) { System.out.println(ex.getMessage()); }
     	}
      return data;
    }
    
    public static ArrayList[] getData(String filename, int ncol) {
    	filename.trim(); String inLine; BufferedReader infile = null; 
     	ArrayList[] data = new ArrayList[ncol]; for (int i=0;i<data.length;i++) { data[i] = new ArrayList(); }
     
     	try {
         	infile = new BufferedReader(new FileReader (filename));  
         	while ((inLine=infile.readLine()) != null) { 	
        		String[] Ldata = inLine.split("\t");
        		for (int i=0;i<ncol;i++) { data[i].add(Ldata[i]); }
         	}
     	} catch (FileNotFoundException ex) { System.out.println("File not found: " + filename); }
     	catch (IOException ex) { System.out.println(ex.getMessage()); }
     	finally { try { if (infile != null) infile.close(); }
         	catch (IOException ex) { System.out.println(ex.getMessage()); }
     	}
      return data;
    }
    
    public void ArrayListOut(ArrayList general) {
    	PrintWriter pw = null; File output = new File(outFileName);
    	try {
        	pw = new PrintWriter(new FileOutputStream(output), true);   	
        	for (int i=0;i<general.size();i++) {
        		pw.print(general.get(i));
      		}
    	}catch (IOException ex) { System.out.println(ex.getMessage()); }
    	finally { if (pw != null) pw.close(); }
    }

    public void Out(ArrayList[] results, boolean ap) {    
    	PrintWriter pw = null;File output = new File(outFileName);    
   		try {
        	pw = new PrintWriter(new FileOutputStream(output, ap), true);
        	for (int i =0;i<results[0].size();i++) {
        		for(int j=0;j<results.length-1;j++) {
        			//System.out.print(results[j].get(i) + "\t");
        			pw.print(results[j].get(i) + "\t"); 
        		}pw.print(results[results.length-1].get(i)+ "\n");
        		//System.out.print(results[results.length-1].get(i)+ "\n");
        	}
    	} catch (IOException ex) {
        	System.out.println(ex.getMessage());
    	} finally { if (pw != null) pw.close(); }
	}

	public static String MakeTable(String filename, double[][] over, String retString) {
        double tot = over.length; double fal = 0; 
        for (int i=0; i<over.length; i++) { if (over[i][0] < 1) fal++; }      
        double tru = tot - fal; double rate = fal / tot * 100;
        retString += filename + "\t" + tot + "\t" + tru + "\t" + fal + "\t" + rate + "\n";     
        return retString;
    }

	public static void MakeFinalReport(String filename, double[][] over, ArrayList[] result) {       
        PrintWriter pw = null; File output = new File(filename);
    	try {
        	pw = new PrintWriter(new FileOutputStream(output), true);   	
        	for (int i=0;i<over.length;i++) {
        		int chr = Integer.parseInt("" + result[0].get(i));
				int start = Integer.parseInt("" + result[1].get(i));
				int stop = Integer.parseInt("" + result[2].get(i));
        		double men = Double.parseDouble("" + result[3].get(i));
        		int num = Integer.parseInt("" + result[4].get(i));
        		double num_over = over[i][0];
        		double per_over = over[i][1];
        		int len = stop-start;		
				if(stop-start > 0) { pw.print(chr + "\t" + start + "\t" + stop + "\t" + men + "\t" + num + "\t" + len + "\t" + num_over + "\t" + per_over  + "\n"); }
      		}
    	}catch (IOException ex) { System.out.println(ex.getMessage()); }
    	finally { if (pw != null) pw.close(); }
    }
	
	public static void Out(int[] chrs, int[] starts, int[] stops, double[] ordata, double[] fildata1, double[] men, String filename) {   
    	PrintWriter pw = null; File output = new File(filename);    
    	try {
        	pw = new PrintWriter(new FileOutputStream(output), true);
        	for (int i =0;i<ordata.length;i++) {
       			pw.print(chrs[i] + "\t" + starts[i] + "\t"+ stops[i] + "\t" + ordata[i] + "\t" + fildata1[i] + "\t" +men[i] + "\n");
        	}     
    	}
    	catch (IOException ex) { System.out.println(ex.getMessage());}
    	finally { if (pw != null) pw.close(); }
	}
	
	public void Out_GFF(ArrayList[] ordata,ArrayList[] calledData, String filename) {    
    	PrintWriter pw = null; File output = new File(filename);
    	try {
        	pw = new PrintWriter(new FileOutputStream(output), true);
        	for (int i =0;i<ordata[0].size();i++) {
        	pw.print("chr" + ordata[0].get(i) + "\t" + "NA" + "\t"+ filename+ "\t" + ordata[1].get(i) + "\t" + ordata[2].get(i) + "\t" + ordata[3].get(i) + "\t" + "." + "\t" + "." + "\t" + "; color 000000" + "\n");
        	}     
        	for (int i =0;i<calledData[0].size();i++) {
        	pw.print("chr" + calledData[0].get(i) + "\t" + "NA" + "\t"+ filename+ "\t" + calledData[1].get(i) + "\t" + calledData[2].get(i) + "\t" + calledData[3].get(i) + "\t" + "." + "\t" + "." + "\t" + "; color FF0000" + "\n");
			}
        
    	} catch (IOException ex) { System.out.println(ex.getMessage()); }
    	finally { if (pw != null) pw.close(); }
	}

}