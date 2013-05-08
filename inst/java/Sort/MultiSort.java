/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package Sort;

import java.util.Arrays;
import java.util.Comparator;
/**
 *
 * @author tf2
 */
public class MultiSort {

    public static String[][] MSort(String[][] toSort, int column) {
        
        Arrays.sort(toSort, new ArrayColumnComparator(column)); 
        return(toSort);
    }

}
class ArrayColumnComparator implements Comparator {

private int columnToSortOn = 0;

    ArrayColumnComparator(int columnToSortOn) {
    this.columnToSortOn = columnToSortOn;
    }

    public int compare(Object o1, Object o2) {

    String[] row1 = (String[])o1;
    String[] row2 = (String[])o2;
    int ret = 0;

    double nb1 = Double.parseDouble(row1[columnToSortOn]);
    double nb2 = Double.parseDouble(row2[columnToSortOn]);


        if ( nb1 == nb2 ) {
            return 0;

        } else if ( nb1 > nb2 ) {
            return 1;

        } else {
            return -1;
        }

   
    }
}

