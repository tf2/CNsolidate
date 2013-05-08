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

public class Clean {

	public static int Min_probe = 5;
	public static double mean_r =0.2;
	public static double Prob_cut =0.2;
	public static int Probe_tol = 3;
	public static String tType = "GADA";
	public static String Gold_standard = "../data/42Mcalls_all_feb09.txt";

	/* Default constructor */
	public Clean() {

	}
	
	public static int[] CleanProbs(double[] data) {	
		int[] passed = new int[data.length]; int count = 0;	
		for(int i=0;i<passed.length;i++) {
			if (Math.abs(data[i]) >Prob_cut) { count = 0; passed[i] = 1;
			} else {	
				count++;
				if(count<Probe_tol) passed[i] = 1;
				else passed[i] = 0;	
			}		
		}	
		return(passed);
	}
		
	public static ArrayList[] MakeInterval(int[] d1, int[] d2, int[] d3, double[] d4, int[] d7) {
		int count = 0; int row_count =0;
		ArrayList[] data = new ArrayList[5];
		for (int i=0;i<data.length;i++) data[i] = new ArrayList();		
		for (int i=0;i<d1.length;i++) {	
			if(d7[i]==1) {	
				count++;	
			} else {
				if (count >0) {
					int staPin = i-count; int stoPin = 0; stoPin = i-1; 
					if (stoPin <0) stoPin=0; 
					if (count > Probe_tol) stoPin= stoPin-Probe_tol; 
					if (Math.abs(d4[staPin]) < mean_r) staPin++;
					if (Math.abs(d4[stoPin]) < mean_r) stoPin--;		
					int ind = stoPin-staPin;
					if (ind < 0) ind = 0;
					double[] v1 = new double[ind]; System.arraycopy(d4,staPin , v1, 0, v1.length);
					double men = 0;
					if (ind > 3) men = exMath.median(v1);
					int len = v1.length; String check = ""+men;		
					if (!check.equals("NaN") &&  !check.equals("NA")) {				
						if (Math.abs(men) >= mean_r && len >= Min_probe) {
							data[0].add(d1[staPin]); data[1].add(d2[staPin]);
							data[2].add(d3[stoPin]);data[3].add(men);
							data[4].add(len); row_count++;		
						}
					}
				}
			count=0;
			}
		}
		return(data);
	}

	public static void ExeShellCommand(String command) {		
        try {		
			System.out.println(command);		
			ProcessBuilder pb = new ProcessBuilder("bash", "-c", command);	
			Process shell = pb.start(); InputStream shellIn = shell.getInputStream();	
            try {
                int shellExitStatus = shell.waitFor();
            } catch (InterruptedException ex) { System.out.println("Problem with shell - waitFor()"); }		
			int c; while ((c = shellIn.read()) != -1) {System.out.write(c);}		
			try {shellIn.close();} catch (IOException ignoreMe) {}	
        } catch(IOException ex) {}
    }
		
	public static void Run(String filename) {
		
		int ind = filename.indexOf(".")-5; String name = filename.substring(0, ind); //String repname = name + "_FinalReport_" + tType + ".txt";
		Handles hand = new Handles(filename, filename); ArrayList[] gDat = hand.getData(6); double[] probs = new double[gDat[0].size()];
		for (int i=0;i<probs.length;i++) probs[i] = Double.parseDouble("" + gDat[5].get(i));
		int[] passed = CleanProbs(probs); 
		int[] d1 = new int[passed.length]; int[] d2 = new int[passed.length]; int[] d3 = new int[passed.length]; double[] d4 = new double[passed.length];
		
		for(int i=0;i<probs.length;i++) {
			d1[i] = Integer.parseInt("" + gDat[0].get(i)); d2[i] = Integer.parseInt("" + gDat[1].get(i));
			d3[i] = Integer.parseInt("" + gDat[2].get(i)); d4[i] = Double.parseDouble("" + gDat[3].get(i));
		}
		
		ArrayList[] result = MakeInterval(d1,d2,d3,d4,passed);			
		int[] chr1 = new int[result[0].size()]; int[] start1 = new int[result[0].size()]; int[] stop1 = new int[result[0].size()];		
		for (int i=0;i<chr1.length;i++) {
			chr1[i] = Integer.parseInt("" + result[0].get(i)); start1[i] = Integer.parseInt("" + result[1].get(i)); stop1[i] = Integer.parseInt("" + result[2].get(i));
		}

		ArrayList[] Gold = Handles.getData(Gold_standard, 3);	
		int[] chr2 = new int[Gold[0].size()]; int[] start2 = new int[Gold[0].size()]; int[] stop2 = new int[Gold[0].size()];	
		for (int i=0;i<chr2.length;i++) {
			chr2[i] = Integer.parseInt("" + Gold[0].get(i)); start2[i] = Integer.parseInt("" + Gold[1].get(i)); stop2[i] = Integer.parseInt("" + Gold[2].get(i));
		}	
		double[][] over = Overlap.ProcessAll(chr1, chr2, start1, start2, stop1, stop2); Handles.MakeFinalReport(filename, over, result);		
		//ExeShellCommand("rm " + filename);
	}
	
	public static void main(String[] args) {	
		String file="";	
		for (int u=0;u<args.length; u++) {		
			int pin = u+1; String s = args[u]; 		
			if (s.equals("-f")) file = args[pin];	
			if (s.equals("-mp")) Min_probe = Integer.parseInt(args[pin]);	
			if (s.equals("-mr")) mean_r = Double.parseDouble(args[pin]);		
			if (s.equals("-pc")) Prob_cut = Double.parseDouble(args[pin]);					
			if (s.equals("-pt")) Probe_tol = Integer.parseInt(args[pin]);		
			if (s.equals("-gold")) Gold_standard = args[pin];		
			if (s.equals("-ty")) tType = args[pin];			
		}
		Run(file);
	}

}

