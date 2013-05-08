import java.io.*;
import java.util.*;
import java.lang.Math.*;
import java.util.ArrayList;
import Sort.*;

public class Gap {

	private static DATAObject DATA;
	public String infile; public String outfile; public Handles hand = null;
	public static String gapFile = "gaps_hg19.txt";

	public Gap(String file, String out) {
		this.infile=file; this.outfile=out;
		this.hand = new Handles(infile, outfile);
	}
 
	public Object[][] runOverGap(String gold, ArrayList[] result) {	
		int[] chr1 = new int[result[0].size()]; int[] start1 = new int[result[0].size()]; int[] stop1 = new int[result[0].size()];
		for (int i=0;i<chr1.length;i++) {
			chr1[i] = Integer.parseInt("" + result[0].get(i)); start1[i] = Integer.parseInt("" + result[1].get(i)); stop1[i] = Integer.parseInt("" + result[2].get(i));
		}
		ArrayList[] Gold = hand.getData(gold, 3);		
		int[] chr2 = new int[Gold[0].size()]; int[] start2 = new int[Gold[0].size()]; int[] stop2 = new int[Gold[0].size()];
		for (int i=0;i<chr2.length;i++) {
			chr2[i] = Integer.parseInt("" + Gold[0].get(i)); start2[i] = Integer.parseInt("" + Gold[1].get(i)); stop2[i] = Integer.parseInt("" + Gold[2].get(i));
		}		
		return( ProcessAllGap(chr1, chr2, start1, start2, stop1, stop2) );
	}

	public int OverOverLap(int chr1, int chr2, int start1, int start2, int stop1, int stop2) {
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

	public Object[][] ProcessAllGap(int[] chr1, int[] chr2, int[] start1, int[] start2, int[] stop1, int[] stop2) {  
		Object[][] over = new Object[chr1.length][3];          
            for (int i=0;i<chr1.length ; i++) {
                int pin = 0; Vector star = new Vector(); Vector stor = new Vector();           	
            	for (int k =0 ;k<chr2.length ; k++) {
                 	int pin2 = OverOverLap(chr1[i], chr2[k], start1[i], start2[k], stop1[i], stop2[k]);
                 	if (pin2 == 1) { 
                 		pin += pin2; star.add(start2[k]); stor.add(stop2[k]); 
                 	} 
             	}
            	if (pin > 0) { 
            		over[i][0] = pin; over[i][1] = Collections.min(star); over[i][2] = Collections.max(stor); 
            		
            	}
                else if (pin == 0) { over[i][0] = 0; over[i][1] = 0; over[i][2] = 0; }
            } 
        return over;
	}

	public void ExeGapFinder() {
		
		this.DATA=new DATAObject(this.infile); Map<String, ArrayList[]> DDATA = this.DATA.DATA; 
		ArrayList[] rep = hand.getData(outfile, 7);	 ArrayList[] newRep = new ArrayList[7];
		for (int i=0;i<newRep.length;i++) newRep[i] = new ArrayList(); 
		
			Object[][] over = runOverGap(gapFile, rep);
			for(int i=0;i<over.length;i++) {	
		
				int rchr1 = Integer.parseInt("" + rep[0].get(i)); int rstart1 = Integer.parseInt("" + rep[1].get(i)); int rstop1 = Integer.parseInt("" + rep[2].get(i)); //Integer.parseInt("" +over[i][1]);			
				int Id1 = Integer.parseInt("" + rep[5].get(i)); int Id2 = Integer.parseInt("" + rep[6].get(i));
		
			if(Integer.parseInt("" + over[i][0]) >0 ) {
				rstart1 = Integer.parseInt("" +over[i][2]); Vector star = new Vector(); Vector rats = new Vector();	Vector rats2 = new Vector();	
				int count=0; int st1=0; int st2=0; 
				for(int j=Id1;j<Id2;j++) {					
					int tchr = Integer.parseInt("" + rep[0].get(i)); int tstart = Integer.parseInt("" + DDATA.get(rep[0].get(i))[0].get(j)); int tstop = Integer.parseInt("" + DDATA.get(rep[0].get(i))[1].get(j));
						if(tstart < Integer.parseInt("" +over[i][1])) {
							rats.add(DDATA.get(rep[0].get(i))[2].get(j)); st1=tstop; count++;				
						} else {
							star.add(tstart);  rats2.add(DDATA.get(rep[0].get(i))[2].get(j));
						}
				}
				
				int p1 = (Id2-(Id2-Id1))+count; int p2 = (Id1+count)+1;
				Object[] aves = new Object[rats.size()];
				for(int k=0;k<aves.length;k++) {
					aves[k] = rats.get(k);
				} double r = exMath.mean(aves);
				
				newRep[0].add(rchr1); newRep[1].add(rep[1].get(i)); newRep[2].add(st1); newRep[3].add(r); newRep[4].add(rats.size()); newRep[5].add(Id1); newRep[6].add(p1);
				rstop1 = Integer.parseInt("" + rep[2].get(i));	 	
				aves = new Object[rats2.size()];
				for(int k=0;k<aves.length;k++) { 
					aves[k] = rats2.get(k); 
				} r = exMath.mean(aves);

				newRep[0].add(rchr1); newRep[1].add(DDATA.get(rep[0].get(i))[0].get(p2)); newRep[2].add(rstop1); newRep[3].add(r); newRep[4].add(rats2.size()); newRep[5].add(p2); newRep[6].add(Id2);
			} else {
				rstop1 = Integer.parseInt("" + rep[2].get(i));
				newRep[0].add(rep[0].get(i)); newRep[1].add(rep[1].get(i)); newRep[2].add(rep[2].get(i)); newRep[3].add(rep[3].get(i)); newRep[4].add(rep[4].get(i)); newRep[5].add(rep[5].get(i)); newRep[6].add(rep[6].get(i));
			}
		}
		
		hand.Out(newRep, false);
	}

 	public static void main(String[] args) {	
		String file =""; String report ="";
		for (int u=0;u<args.length; u++) {	
			int pin = u+1; String s = args[u]; 
			if (s.equals("-f")) file = args[pin];
			if (s.equals("-r")) report = args[pin];
			if (s.equals("-g")) gapFile = args[pin];
		}
		Gap gap = new Gap(file, report); gap.ExeGapFinder();
	}

}





