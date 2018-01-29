/*----------------------------------------------------------------
 *  Author:        Junjiajia Long
 *  Written:       3/14/2015
 *  Last updated:  3/16/2015
 *
 *  Compilation:   javac KdTree.java
 *  Execution:     java KdTree
 *----------------------------------------------------------------*/
import java.util.Iterator;

public class KdTree {
    private Node root; // root node of the BST
    private int N; // total number of nodes in the BST
    
    public KdTree() {                               // construct an empty tree of nodes
        root = null;
        N = 0;
    }
    public boolean isEmpty() {                        // is the set empty?
        return N == 0;
    }
    public int size() {                               // number of points in the set
        return N;
    }
    
    public void insert(Point2D p) {                   // add the point to the set (if it is not already in the set)
        if (p == null) throw new NullPointerException("insert() argument null");
        if (root == null) { // create the first node
            root = new Node(p, true, new RectHV(0.0, 0.0, 1.0, 1.0));
            N = 1;
            return;
        }
        Node curr = root;
        Node last = null;
        double compare = 0; // result of comparison
        while (curr != null) { // find the place to insert
            if (curr.xSep) compare = p.x() - curr.point.x();
            else compare = p.y() - curr.point.y();
            last = curr;
            if (compare < 0) curr = curr.left;
            else if (compare > 0) curr = curr.right;
            else if (!p.equals(curr.point)) curr = curr.right;
            else return; // exit when the point is found as duplicate
        }
        
        // find the box containing the point and create link
        RectHV subBox;
        if (compare < 0) {
            if (last.xSep) subBox = new RectHV(last.box.xmin(), last.box.ymin(), last.point.x(), last.box.ymax());
            else subBox = new RectHV(last.box.xmin(), last.box.ymin(), last.box.xmax(), last.point.y());
            curr = new Node(p, !last.xSep, subBox);
            last.left = curr;
        } else {
            if (last.xSep) subBox = new RectHV(last.point.x(), last.box.ymin(), last.box.xmax(), last.box.ymax());
            else subBox = new RectHV(last.box.xmin(), last.point.y(), last.box.xmax(), last.box.ymax());
            curr = new Node(p, !last.xSep, subBox);
            last.right = curr;
        }
        N = N + 1;
    }
    public boolean contains(Point2D p) {              // does the set contain point p?
        if (p == null) throw new NullPointerException("contains() argument null");
        Node curr = root;
        while (curr != null) {
            double compare; // result of comparison
            if (curr.xSep) compare = p.x() - curr.point.x();
            else compare = p.y() - curr.point.y();
            if (compare < 0) curr = curr.left;
            else if (compare > 0) curr = curr.right;
            else if (!p.equals(curr.point)) curr = curr.right;
            else return true; // the point is found
        }
        return false;
    }
    
    public void draw() {                              // draw all points to standard draw
        if (root != null) drawAllChildren(root);
    }
    private void drawAllChildren(Node curr) { // draw the node point with splitting line
        // draw the point
        StdDraw.setPenColor(StdDraw.BLACK);
        StdDraw.setPenRadius(.01);
        curr.point.draw();
        
        // draw the splitting line defined by box and defint the subboxes
        StdDraw.setPenRadius();
        if (curr.xSep) {
            StdDraw.setPenColor(StdDraw.RED);
            StdDraw.line(curr.point.x(), curr.box.ymin(), curr.point.x(), curr.box.ymax());
        } else {
            StdDraw.setPenColor(StdDraw.BLUE);
            StdDraw.line(curr.box.xmin(), curr.point.y(), curr.box.xmax(), curr.point.y());
        }
        
        // recurse to subtrees
        if (curr.left != null) drawAllChildren(curr.left);
        if (curr.right != null) drawAllChildren(curr.right);
    }
    
    public Iterable<Point2D> range(RectHV rect) {     // all points that are inside the rectangle
        if (rect == null) throw new NullPointerException("range() argument null");
        Stack<Point2D> insiders = new Stack<Point2D>();
        if (root != null) insiders = rQuery(rect, root, insiders);
        return insiders;
    }
    private Stack<Point2D> rQuery(RectHV rct, Node curr, Stack<Point2D> in) { // check the current node and query the subtrees
        if (rct.contains(curr.point)) in.push(curr.point);
        if (curr.left != null && rct.intersects(curr.left.box)) in = rQuery(rct, curr.left, in);
        if (curr.right != null && rct.intersects(curr.right.box)) in = rQuery(rct, curr.right, in);
        return in;
    }
    
    public Point2D nearest(Point2D p) {               // a nearest neighbor in the set to point p, null if the set is empty
        if (p == null) throw new NullPointerException("nearest() argument null");
        double nearestDistanceSquared = Double.POSITIVE_INFINITY;
        Point2D nearestPoint = null;
        if (root != null) nearestPoint = nQuery(p, root, root.point);
        return nearestPoint;
    }
    private Point2D nQuery(Point2D p, Node curr, Point2D nP) { // check the current node and query the subtrees
        double cDS = p.distanceSquaredTo(curr.point);
        double nDS = p.distanceSquaredTo(nP);
        if (cDS < nDS) nP = curr.point;
        // a node is searched only if it might contain a point that is closer than the best one found so far
        boolean leftSubtreeExplore = curr.left != null && curr.left.box.distanceSquaredTo(p) < nDS;
        boolean rightSubtreeExplore = curr.right != null && curr.right.box.distanceSquaredTo(p) < nDS;
        if (leftSubtreeExplore && rightSubtreeExplore) {
            // determine priority of subtree exploration
            boolean leftPriority;
            if (curr.xSep) {
                if (p.x() < curr.point.x()) leftPriority = true;
                else leftPriority = false;
            } else {
                if (p.y() < curr.point.y()) leftPriority = true;
                else leftPriority = false;
            }
            // explore according to priority
            if (leftPriority) {
                nP = nQuery(p, curr.left, nP);
                nP = nQuery(p, curr.right, nP);
            } else {
                nP = nQuery(p, curr.right, nP);
                nP = nQuery(p, curr.left, nP);
            }
        } else if (leftSubtreeExplore) {
            nP = nQuery(p, curr.left, nP);
        } else if (rightSubtreeExplore) {
            nP = nQuery(p, curr.right, nP);
        }
        return nP;
    }

    private class Node { // Node in the BST
        public final Point2D point;
        public Node left;
        public Node right;
        public final boolean xSep;
        public final RectHV box;
        public Node(Point2D p, boolean sep, RectHV b) {
            point = p;
            left = null;
            right = null;
            xSep = sep;
            box = b;
        }
    }
    
    public static void main(String[] args) {          // unit testing of the methods (optional) 
        String filename = args[0];
        In in = new In(filename);

        // initialize the two data structures with point from standard input
        KdTree kdtree = new KdTree();
        while (!in.isEmpty()) {
            double x = in.readDouble();
            double y = in.readDouble();
            Point2D p = new Point2D(x, y);
            kdtree.insert(p);
        }
        
        // draw all points
        kdtree.draw();
    }
}