/*----------------------------------------------------------------
 *  Author:        Junjiajia Long
 *  Written:       2/21/2015
 *  Last updated:  2/24/2015
 *
 *  Compilation:   javac Brute.java
 *  Execution:     java Brute
 *----------------------------------------------------------------*/
import java.util.Arrays;

public class Brute {
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
        
        // loop over all 4-tuples
        for (int i = 0; i < N-3; i++) {
            for (int j = i+1; j < N-2; j++) {
                for (int k = j+1; k < N-1; k++) {
                    if (p[i].slopeTo(p[j]) == p[i].slopeTo(p[k])) {
                        for (int l = k+1; l < N; l++) {
                            if (p[i].slopeTo(p[j]) == p[i].slopeTo(p[l])) {
                                StdOut.print(p[i].toString());
                                StdOut.print(" -> ");
                                StdOut.print(p[j].toString());
                                StdOut.print(" -> ");
                                StdOut.print(p[k].toString());
                                StdOut.print(" -> ");
                                StdOut.println(p[l].toString());
                                p[i].drawTo(p[l]);
                            }
                        }
                    }
                }
            }
        }
    }
}