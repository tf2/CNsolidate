/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package Sort;

/**
 *
 * @author tf2
 */

import java.util.Arrays;

/**
 * @author Ramakrishnan Sundararaj
 *
 */
public class SortManager
{
	private static int currCol = 0;
	public static final int ASC = 1;
	public static final int DES = -1;
	private static int currOrd = ASC;
	private static int length = 0;
	private static int width = 0;

	public static int getCurrCol()
	{
		return currCol;
	}

	public static String[][] sort(String[][] twodime, int[] column, int[] order)
	{
		length = twodime.length;
		width = length > 0? twodime[0].length: 0;

		Object[] arr = new Object[length];
		for(int i=0; i<length; i++)
		{
			arr[i] = new SortElem(twodime[i]);
		}

		for(int i=0, len=column.length; i<len; i++)
		{
			currCol = column[i];
			currOrd = order[i];
			Arrays.sort(arr);
			for(int idx=arr.length-1; idx >= 0; idx-- )
			{
				((SortElem) arr[idx]).setNext(currCol);
			}
		}

		String[][] retArr = getTwoDimeArr(arr);

		return retArr;
	}

	public static int compare(String[] row1, String[] row2) {
		
        int ret = 0;
        double nb1 = Double.parseDouble(row1[currCol]);
        double nb2 = Double.parseDouble(row2[currCol]);

            if ( nb1 == nb2 ) {
                return 0;
            } else if ( nb1 > nb2 ) {
                return 1;
            } else {
                return -1;
            }

	}

	public static String[][] getTwoDimeArr(Object[] arr)
	{
		String[][] retArr = new String[length][width];
		for(int i=0; i<length; i++) {
			SortElem elem = (SortElem) arr[i];
			retArr[i] = elem.getRow();
		}
		return retArr;
	}

	public static void main(String args[])
	{

            //double[] t = Handles.Handles.InRatNimblegen("/Users/tf2/Desktop/369909A03_S01_segMNT.txt", Handles.Handles.ArraySizer("/Users/tf2/Desktop/369909A03_S01_segMNT.txt"));
           // String[][] myS = waveHandle.WaveNimblegen("/Users/tf2/Desktop/369909A03_S01_segMNT.txt", t, Handles.Handles.ArraySizer("/Users/tf2/Desktop/369909A03_S01_segMNT.txt"));
           // System.out.println("read");
          //  myS = sort(myS, new int[]{0, 1}, new int[]{ASC, ASC}); // sort(2-dime-array, columns, order);
            
            //System.out.println("sorted");
            
          //  double[] t = Handles.Handles.InAll("/Users/tf2/Desktop/US22502573_252347210003_S01_CGH_105_Dec08.txt", Handles.Handles.ArraySizer("/Users/tf2/Desktop/US22502573_252347210003_S01_CGH_105_Dec08.txt"));
         //   double[] tt = waveHandle.WaveAgilent("/Users/tf2/Desktop/US22502573_252347210003_S01_CGH_105_Dec08.txt", t, Handles.Handles.ArraySizer("/Users/tf2/Desktop/US22502573_252347210003_S01_CGH_105_Dec08.txt"));
           
          //  for (int i=0;i<tt.length;i++) {
         //       System.out.println(t[i] + "\t" + tt[i]);
          //  }


            /*System.out.println("read");
            myS = sort(myS, new int[]{0, 1}, new int[]{ASC, ASC}); // sort(2-dime-array, columns, order);
            
            System.out.println("sorted");
            
            
            
            for(int i=0, len=myS.length; i<len; i++)
		{
			String[] arr = myS[i];
			for(int j=0, wid=arr.length; j<wid; j++) {
				System.out.print(arr[j] + " | ");
			}
			System.out.println("\n------------------------------------------------------");
		}
            for ( int i=0;i<100;i++) {
            System.out.println("\n------------------------------------------------------");
            }
            
            
             myS = sort(myS, new int[]{2}, new int[]{ASC});
            
             
             for(int i=0, len=myS.length; i<len; i++)
		{
			String[] arr = myS[i];
			for(int j=0, wid=arr.length; j<wid; j++) {
				System.out.print(arr[j] + " | ");
			}
			System.out.println("\n------------------------------------------------------");
		}
             
             
            /*String[][] twodime = new String[][]{{"1","20","5"},
				{"5","11","6"},
				{"1","3","7"},
				{"2","25","8"},
				{"6","22","8"},
				{"4","12","21"}
				};


               
		twodime = sort(twodime, new int[]{0, 1}, new int[]{ASC, ASC}); // sort(2-dime-array, columns, order);

		for(int i=0, len=twodime.length; i<len; i++)
		{
			String[] arr = twodime[i];
			for(int j=0, wid=arr.length; j<wid; j++) {
				System.out.print(arr[j] + " | ");
			}
			System.out.println("\n------------------------------------------------------");
		}
               */
	}
}

