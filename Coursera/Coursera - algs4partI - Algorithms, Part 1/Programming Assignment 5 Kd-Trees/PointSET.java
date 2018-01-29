/*----------------------------------------------------------------
 *  Author:        Junjiajia Long
 *  Written:       3/9/2015
 *  Last updated:  3/15/2015
 *
 *  Compilation:   javac PointSET.java
 *  Execution:     java PointSET
 *----------------------------------------------------------------*/
import java.util.Iterator;

public class PointSET {
    private SET<Point2D> points; // set of all points
    
    public PointSET() {                               // construct an empty set of points
        points = new SET<Point2D>();
    }
    public boolean isEmpty() {                        // is the set empty?
        return points.isEmpty();
    }
    public int size() {                               // number of points in the set
        return points.size();
    }
    public void insert(Point2D p) {                   // add the point to the set (if it is not already in the set)
        if (p == null) throw new NullPointerException("insert() argument null");
        points.add(p);
    }
    public boolean contains(Point2D p) {              // does the set contain point p?
        if (p == null) throw new NullPointerException("contains() argument null");
        return points.contains(p);
    }
    public void draw() {                              // draw all points to standard draw
        Iterator<Point2D> allPoints = points.iterator();
        while (allPoints.hasNext()) allPoints.next().draw();
    }
    public Iterable<Point2D> range(RectHV rect) {     // all points that are inside the rectangle
        if (rect == null) throw new NullPointerException("range() argument null");
        Stack<Point2D> insiders = new Stack<Point2D>();
        Iterator<Point2D> allPoints = points.iterator();
        while (allPoints.hasNext()) {
            Point2D in = allPoints.next();
            if (rect.contains(in)) insiders.push(in);
        }
        return insiders;
    }
    public Point2D nearest(Point2D p) {               // a nearest neighbor in the set to point p; null if the set is empty
        if (p == null) throw new NullPointerException("nearest() argument null");
        double nearestDistanceSquared = 2;
        Point2D nearestPoint = null;
        Iterator<Point2D> allPoints = points.iterator();
        while (allPoints.hasNext()) {
            Point2D nb = allPoints.next();
            double currentDistanceSquared = p.distanceSquaredTo(nb);
            if ( currentDistanceSquared < nearestDistanceSquared) {
                nearestDistanceSquared = currentDistanceSquared;
                nearestPoint = nb;
            }
        }
        return nearestPoint;
    }

    public static void main(String[] args) {          // unit testing of the methods (optional) 
    }
}