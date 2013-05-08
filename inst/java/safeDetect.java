import java.io.*;
import java.util.*;
import java.lang.Math.*;
import java.util.ArrayList;


public class safeDetect {

	private double minRat; private int Pnum;
	private String datafile=""; private String reportfile="";

	public safeDetect(String datafile, double minRat, int Pnum) {
		this.datafile = datafile; this.minRat=minRat; this.Pnum=Pnum;
		int index=datafile.lastIndexOf("."); String inName=datafile; reportfile=datafile.substring(0, index)+"_FinalReport_safe.txt";	
		safelyDetect(datafile, reportfile);			
	}

    public void safelyDetect(String inname, String outname) {    
    	Handles hand = new Handles(inname, outname); ArrayList[] data = hand.getData(4);
    	ArrayList res = new ArrayList();
 		Object[] chrs = new Object[data[0].size()];
 		for(int i=0;i<chrs.length;i++) chrs[i] = data[0].get(i);		
 		ArrayList uni = new ArrayList(); java.util.List lis = Arrays.asList(chrs);
		for (Object x : lis) { if (!uni.contains(x)) uni.add(x); }   		
   		ArrayList[] ratio = new ArrayList[uni.size()]; ArrayList[] starts = new ArrayList[uni.size()]; ArrayList[] stops = new ArrayList[uni.size()];		
   		for(int i=0;i<ratio.length;i++) { ratio[i] = new ArrayList(); starts[i] = new ArrayList(); stops[i] = new ArrayList(); }
   		
   		int index = 0; ratio[index].add(data[3].get(0)); starts[index].add(data[1].get(0)); stops[index].add(data[2].get(0)); 		
   		for(int i=1;i<data[0].size();i++) {
  			int loc1 = Integer.parseInt(""+data[0].get(i)); int loc2 = Integer.parseInt(""+data[0].get(i-1));
  			if(loc1 != loc2) index++;
   			ratio[index].add(data[3].get(i)); starts[index].add(data[1].get(i)); stops[index].add(data[2].get(i));
   		}		
   		
   		for (int i=0;i<uni.size();i++) {	   		
   		   	int check1 = 0; int check2 = 0;	   	
   			for(int j=0;j<ratio[i].size()-Pnum;j++) {		
   				int k = j; int kk = j; Object[] rat = ratio[i].toArray();
				
   				while(Double.parseDouble("" + ratio[i].get(k)) > minRat && k<ratio[i].size()-1) { check1++; k++; }		
   				while(Double.parseDouble("" + ratio[i].get(kk)) < -minRat && kk<ratio[i].size()-1) { check2++; kk++; }   
   			
   				if (check1 > Pnum) {
   					int ch = Integer.parseInt(""+uni.get(i)); Object[] ar = new Object[check1-1];
   					System.arraycopy(rat, j , ar, 0, ar.length); double men = exMath.mean(ar);
					if (men > minRat) res.add(ch + "\t" + starts[i].get(j)  + "\t" + stops[i].get(j+check1-1) + "\t" + men + "\t" + ar.length + "\n");
   				 	j = j+check1;
   				}			
   				if (check2 > Pnum) {
   					int ch = Integer.parseInt(""+uni.get(i)); Object[] ar = new Object[check2];
   					System.arraycopy(rat, j , ar, 0, ar.length); double men = exMath.mean(ar);   					
   					if (men < -minRat) res.add(ch + "\t" + starts[i].get(j)  + "\t" + stops[i].get(j+check2-1) + "\t" + men + "\t" + ar.length + "\n");
   					j = j+check2;
   				}			
   				check1 = 0; check2 = 0;
   			}
   		} hand.ArrayListOut(res);
    }
    
   public static void main(String[] args) {	
		String file=""; double mr=0.3; int pn=3;		
		for (int u=0;u<args.length; u++) {			
			int pin = u+1;
            String s = args[u]; 			
			if (s.equals("-f")) file = args[pin];
			if(s.equals("-mr")) mr = Double.parseDouble("" + args[pin]);
			if(s.equals("-np")) pn = Integer.parseInt("" + args[pin]);
		} safeDetect sd = new safeDetect(file, mr, pn);
	}
    
}






