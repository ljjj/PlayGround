/*----------------------------------------------------------------
 *  Author:        Junjiajia Long
 *  Written:       2/19/2015
 *  Last updated:  2/21/2015
 *
 *  Compilation:   javac Point.java
 *  Execution:     java Point
 *  Dependencies:  StdDraw.java
 *----------------------------------------------------------------*/
import java.util.Arrays;
import java.util.Comparator;

public class Point implements Comparable<Point> {
    public final Comparator<Point> SLOPE_ORDER = new SlopeOrder(); // compare points by slope to this point

    private final int x;                              // x coordinate
    private final int y;                              // y coordinate
    
    public Point(int x, int y) {                       // construct the point (x, y)
        this.x = x;
        this.y = y;
    }
    
    public   void draw() {                             // draw this point
        StdDraw.point(x, y);
    }
    
    public   void drawTo(Point that) {                 // draw the line segment from this point to that point
        StdDraw.line(this.x, this.y, that.x, that.y);
    }
    
    public String toString() {                         // string representation
        return "(" + x + ", " + y + ")";
    }

    public    int compareTo(Point that) {              // is this point lexicographically smaller than that point?
        if (this.y < that.y) { return -1;}
        else if (this.y > that.y) { return 1;}
        else if (this.x < that.x) { return -1;}
        else if (this.x > that.x) { return 1;}
        else { return 0;}
    }
    
    public double slopeTo(Point that) {                // the slope between this point and that point
        if (this.x == that.x) {
            if (this.y == that.y) { return Double.NEGATIVE_INFINITY;}
            else { return Double.POSITIVE_INFINITY;}
        } else if (this.y == that.y) {return (double)+0.0;}
        return ((double) that.y - this.y)/(that.x - this.x);
    }
    
    private class SlopeOrder implements Comparator<Point> {
        public int compare(Point here, Point there) {
            if (slopeTo(here) < slopeTo(there)) {return -1;}
            else if (slopeTo(here) > slopeTo(there)) {return 1;}
            else {return 0;}
        }
    }
    
    // unit test
    public static void main(String[] args) {
        // create 6 points for testing
        Point[] P = new Point[6];
        P[0] = new Point(0,0);
        P[1] = new Point(1,0);
        P[2] = new Point(1,1);
        P[3] = new Point(0,1);
        P[4] = new Point(-1,0);
        P[5] = new Point(0,-1);
        
        // draw points and lines
        StdDraw.setXscale(-2, 2);
        StdDraw.setYscale(-2, 2);
        P[0].draw();
        P[1].draw();
        P[2].draw();
        P[3].draw();
        P[4].draw();
        P[5].draw();
        P[0].drawTo(P[2]);
        P[5].drawTo(P[4]);
        
        // compare points for output
        Arrays.sort(P, 0, 6);
        StdOut.println("Sort by natural order:");
        for (int i = 0; i < 6; i++) {StdOut.println(P[i].toString());}
        
        // check slope
        StdOut.println("Slope to (0, 0):");
        for (int i = 0; i < 6; i++) {
            StdOut.print(P[i].toString());
            StdOut.print(":");
            StdOut.println(P[2].slopeTo(P[i]));}
        
        // check comparator
        Arrays.sort(P, 0, 6, P[2].SLOPE_ORDER);
        StdOut.println("Sort by (0, 0):");
        for (int i = 0; i < 6; i++) {StdOut.println(P[i].toString());}
    }
}