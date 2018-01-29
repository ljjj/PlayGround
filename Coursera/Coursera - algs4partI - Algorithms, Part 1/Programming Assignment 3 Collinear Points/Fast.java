/*----------------------------------------------------------------
 *  Author:        Junjiajia Long
 *  Written:       2/24/2015
 *  Last updated:  2/24/2015
 *
 *  Compilation:   javac Fast.java
 *  Execution:     java Fast
 *----------------------------------------------------------------*/
import java.util.Arrays;

public class Fast {
    public static void main(String[] args) {
        // read in and sort
        In in = new In(args[0]);      // input file
        int N = in.readInt();         // number of points
        Point[] p = new Point[N];
        StdDraw.setXscale(0, 32768);
        StdDraw.setYscale(0, 32768);
        for (int i = 0; i < N; i++) {
            int x = in.readInt();
            int y = in.readInt();
            p[i] = new Point(x,y);
            p[i].draw();
        }
        Arrays.sort(p, 0, N);
        
        for (int i = 0; i < N-3; i++) {
            // for each point, copy a new point array and sort according to slope
            Point[] q = new Point[N];
            for (int c = 0; c < N; c++) { q[c] = p[c]; }
            Arrays.sort(q, 0, N, p[i].SLOPE_ORDER);
            
            int j = 1;
            while (j < N-2) {
                double firstSlope = q[0].slopeTo(q[j]); // first other point to look at slope
                boolean newLine = q[j].compareTo(q[0]) > 0; // check if the line is already found
                int k = j+1;
                while ((k < N) && (q[0].slopeTo(q[k]) == firstSlope)) {
                    newLine = newLine && (q[k].compareTo(q[0]) > 0);
                    k++; // search for the next few points with the same slope
                }
                if ((k > j+2) && newLine) { // output if more than 3 other points are found in a new line
                    StdOut.print(q[0].toString());
                    for (int l = j; l < k; l++) {
                        StdOut.print(" -> ");
                        StdOut.print(q[l].toString());
                    }
                    StdOut.println();
                    q[0].drawTo(q[k-1]);
                }
                j = k; // start with the next point with different slope
            }
        }
    }
}