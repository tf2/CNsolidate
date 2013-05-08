import java.io.*;
import java.util.*;
import java.lang.Math.*;
import java.util.ArrayList;

public class fiFast {
	
	private String fileName = ""; private String outName = "";
	
	public fiFast(String file, String out) {
		this.fileName = file; this.outName=out;
	}

	public void extractFcalls() {		
		Handles hand = new Handles(fileName, outName);
		ArrayList[] data = hand.getData(6);
		ArrayList res = new ArrayList();
		ArrayList[] AA = new ArrayList[2];
		AA[0] = new ArrayList(); AA[1] = new ArrayList();
		Object[] data2 = data[3].toArray();		
		int pin = 0; int p1 = 0; int pp=0; int x=0;
		while( pp <= data[0].size()-1 ) {	
			if(pin >=data[4].size()-1) { pp = data[0].size(); break; }
			int c = Integer.parseInt( "" + data[4].get(x) );
			if(c==2) {
				while(c!=0) {
				if(pin >=data[4].size()-1) { pp = data[0].size(); break; }
					c = Integer.parseInt( "" + data[4].get(pin) ); pin++;		
				}
				AA[0].add(x); AA[1].add(pin-1);
				if(pin-1 > x) {
					Object[] temp = new Object[(pin-1)-x]; System.arraycopy(data2, x, temp, 0, temp.length); double men = exMath.mean(temp);
					// This algorithm is far too sensitive!!!
					if( temp.length>2 & Math.abs(men) >0.3 & Integer.parseInt(""+data[1].get(x)) < Integer.parseInt(""+data[2].get(pin-1)) ) {
						res.add(data[0].get(x) + "\t" + data[1].get(x) + "\t" + data[2].get(pin-1) + "\t" +  men + "\t" + temp.length + "\n");
					}
				} p1++; x = pin;
			} pin=x; x++; pp++;
		} hand.ArrayListOut(res);
	}

	public static void main(String[] args) {
		String file =""; String out ="";	
		for (int u=0;u<args.length; u++) {	
			int pin = u+1; String s = args[u]; 	
			if (s.equals("-f")) file = args[pin];
			if (s.equals("-o")) out = args[pin];
		}
		fiFast fif = new fiFast(file, out); fif.extractFcalls();	
	}

}
