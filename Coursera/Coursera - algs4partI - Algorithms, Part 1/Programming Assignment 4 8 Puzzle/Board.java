/*----------------------------------------------------------------
 *  Author:        Junjiajia Long
 *  Written:       3/7/2015
 *  Last updated:  3/7/2015
 *
 *  Compilation:   javac Board.java
 *  Execution:     java Board puzzle00.txt
 *----------------------------------------------------------------*/

public class Board {
    private final char[] tiles; // block values, 0 represents blank
    private final int N; // board dimension < 128
    private final int ham; // Hamming distance
    private final int man; // Manhattan distance
    private final int blank; // location of blank
    
    public Board(int[][] blocks) {         // construct a board from an N-by-N array of blocks
        N = blocks.length;                 // (where blocks[i][j] = block in row i, column j)
        // initialize tiles
        tiles = new char[N*N];
        int b = -1;
        int hamming = 0;
        int manhattan = 0;
        
        // read tiles and calculate distances
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                int value = blocks[i][j];
                tiles[N*i+j] = (char) value;
                if (value == 0) b = N*i+j;
                else {
                    if (value != N*i+j+1) hamming++; // calculate Hamming distance
                    manhattan += mandis(i,j,value); // calculate Manhattan distance
                }
            }
        }
        blank = b;
        ham = hamming;
        man = manhattan;
    }
    public int dimension() {               // board dimension N
        return N;
    }
    private Board(char[] nbtiles, int nbN, int nbham, int nbman, int nbblank) { // for neighbors and twins
        N = nbN;
        tiles = new char[N*N];
//        int testham = 0;
//        int testman = 0;
//        int testblank = -1;
        for (int k = 0; k < N*N; k++) {
            tiles[k] = nbtiles[k];
//            int value = (int)tiles[k];
//            if (value == 0) testblank = k;
//            else {
//                if (value != k+1) testham++; // calculate Hamming distance
//                testman += mandis(k/N,k%N,value); // calculate Manhattan distance
//            }
        }
//        assert testham == nbham;
//        assert testman == nbman;
//        assert testblank == nbblank;
        ham = nbham;
        man = nbman;
        blank = nbblank;
    }
    
    public int hamming() {                 // number of blocks out of place
        return ham;
    }
    public int manhattan() {               // sum of Manhattan distances between blocks and goal
        return man;
    }
    private int mandis(int i, int j, int value) { // Manhattan distance of a single tile
        int supposedI = (value-1)/N;
        int supposedJ = (value-1)%N;
        return Math.abs(i - supposedI) + Math.abs(j - supposedJ);
    }
    
    public boolean isGoal() {              // is this board the goal board?
        return ham == 0;
    }
    public Board twin() {                  // a board that is obtained by exchanging two adjacent blocks in the same row
        // copy the board
//        int[][] temp = copy();
        // exchange two adjacent blocks
        Board twin;
        int newham = ham;
        int newman = man;
        if (blank >= N) { // if blank is not in the first row, swap first two tiles
//            temp[0][0] = (int)tiles[1];
//            temp[0][1] = (int)tiles[0];
            char temp = tiles[0];
            tiles[0] = tiles[1];
            tiles[1] = temp;
            if((char) tiles[1] == 1) newham++;
            if((char) tiles[0] == 2) newham++;
            if((char) tiles[0] == 1) newham--;
            if((char) tiles[1] == 2) newham--;
            newman -= mandis(0,0,(char) tiles[1]);
            newman -= mandis(0,1,(char) tiles[0]);
            newman += mandis(0,0,(char) tiles[0]);
            newman += mandis(0,1,(char) tiles[1]);
            twin = new Board(tiles, N, newham, newman, blank);
            tiles[1] = tiles[0];
            tiles[0] = temp;
        }
        else { // swap the first two tiles in the second row
//            temp[1][0] = (int)tiles[N+1];
//            temp[1][1] = (int)tiles[N];
            char temp = tiles[N];
            tiles[N] = tiles[N+1];
            tiles[N+1] = temp;
            if((char) tiles[N+1] == N+1) newham++;
            if((char) tiles[N] == N+2) newham++;
            if((char) tiles[N] == N+1) newham--;
            if((char) tiles[N+1] == N+2) newham--;
            newman -= mandis(1,0,(char) tiles[N+1]);
            newman -= mandis(1,1,(char) tiles[N]);
            newman += mandis(1,0,(char) tiles[N]);
            newman += mandis(1,1,(char) tiles[N+1]);
            twin = new Board(tiles, N, newham, newman, blank);
            tiles[N+1] = tiles[N];
            tiles[N] = temp;
        }
        return twin;
    }
    public boolean equals(Object y) {      // does this board equal y?
        if (y == this) return true;
        if (y == null) return false;
        if (y.getClass() != this.getClass()) return false;
        Board that = (Board) y;
        boolean eq = true;
        for (int k = 0; k < N*N; k++)
            eq = eq && tiles[k] == that.tiles[k];
        return eq;
    }
    
    public Iterable<Board> neighbors() {   // all neighboring boards
        Stack<Board> nb = new Stack<Board>(); // stack of neighbors
        // create flags to tell whether a move is valid
        boolean bUable = true; // blank can move up
        boolean bDable = true; // blank can move down
        boolean bLable = true; // blank can move left
        boolean bRable = true; // blank can move right
        
        // set the flags depending on the blank position
        if (blank == 0) { // top-left corner
            bUable = false;
            bLable = false;
        }
        if (blank == N-1) { // top-right corner
            bUable = false;
            bRable = false;
        }
        if (blank == N*N-N) { // bottom-left corner
            bDable = false;
            bLable = false;
        }
        if (blank == N*N-1) { // bottom-right corner
            bDable = false;
            bRable = false;
        }
        if (blank < N) bUable = false; // top row
        if (blank >= N*N-N) bDable = false; // bottom row
        if (blank%N == 0) bLable = false; // leftmost column
        if (blank%N == N-1) bRable = false; // rightmost column
        
        // push neighbors to stack depending on flags
        if (bUable) nb.push(blankUp());
        if (bDable) nb.push(blankDown());
        if (bLable) nb.push(blankLeft());
        if (bRable) nb.push(blankRight());
        return nb;
    }
    
    public String toString() {             // string representation of this board (in the output format specified below)
        StringBuilder s = new StringBuilder();
        s.append(N + "\n");
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                s.append(String.format("%2d ", (int)tiles[N*i+j]));
            }
            s.append("\n");
        }
        return s.toString();
    }
    
//    private int[][] copy() {                 // create an array as copy
//        int[][] cp = new int[N][N];
//        for (int i = 0; i < N; i++)
//            for (int j = 0; j < N; j++)
//                cp[i][j] = (int)tiles[N*i+j];
//        assert cp[blank/N][blank%N] == 0;
//        return cp;
//    }
    private Board blankUp() {              // move the blank up
        // check validity of the move
        if (blank < N)
            throw new IndexOutOfBoundsException
            ("Trying to move blank up in top row");
        // move the blank
//        int[][] temp = copy();
//        int blankI = blank/N;
//        int blankJ = blank%N;
//        temp[blankI][blankJ] = (int)tiles[blank-N];
//        temp[blankI-1][blankJ] = 0;
        int newham = ham;
        int newman = man;
        tiles[blank] = tiles[blank-N];
        tiles[blank-N] = 0;
        if((char) tiles[blank] == blank-N+1) newham++;
        if((char) tiles[blank] == blank+1) newham--;
        newman -= mandis(blank/N-1, blank%N, (char) tiles[blank]);
        newman += mandis(blank/N, blank%N, (char) tiles[blank]);
        Board nb = new Board(tiles, N, newham, newman, blank-N);
        tiles[blank-N] = tiles[blank];
        tiles[blank] = 0;
        return nb;
    }
    private Board blankDown() {              // move the blank down
        // check validity of the move
        if (blank >= N*N-N)
            throw new IndexOutOfBoundsException
            ("Trying to move blank down in bottom row");
        // move the blank
//        int[][] temp = copy();
//        int blankI = blank/N;
//        int blankJ = blank%N;
//        temp[blankI][blankJ] = (int)tiles[blank+N];
//        temp[blankI+1][blankJ] = 0;
//        Board mv = new Board(temp);
//        return mv;
        int newham = ham;
        int newman = man;
        tiles[blank] = tiles[blank+N];
        tiles[blank+N] = 0;
        if((char) tiles[blank] == blank+N+1) newham++;
        if((char) tiles[blank] == blank+1) newham--;
        newman -= mandis(blank/N+1, blank%N, (char) tiles[blank]);
        newman += mandis(blank/N, blank%N, (char) tiles[blank]);
        Board nb = new Board(tiles, N, newham, newman, blank+N);
        tiles[blank+N] = tiles[blank];
        tiles[blank] = 0;
        return nb;
    }
    private Board blankLeft() {              // move the blank left
        // check validity of the move
        if (blank%N == 0)
            throw new IndexOutOfBoundsException
            ("Trying to move blank left in leftmost column");
        // move the blank
//        int[][] temp = copy();
//        int blankI = blank/N;
//        int blankJ = blank%N;
//        temp[blankI][blankJ] = (int)tiles[blank-1];
//        temp[blankI][blankJ-1] = 0;
//        Board mv = new Board(temp);
//        return mv;
        int newham = ham;
        int newman = man;
        tiles[blank] = tiles[blank-1];
        tiles[blank-1] = 0;
        if((char) tiles[blank] == blank) newham++;
        if((char) tiles[blank] == blank+1) newham--;
        newman -= mandis((blank-1)/N, (blank-1)%N, (char) tiles[blank]);
        newman += mandis(blank/N, blank%N, (char) tiles[blank]);
        Board nb = new Board(tiles, N, newham, newman, blank-1);
        tiles[blank-1] = tiles[blank];
        tiles[blank] = 0;
        return nb;
    }
    private Board blankRight() {              // move the blank right
        // check validity of the move
        if (blank%N == N-1)
            throw new IndexOutOfBoundsException
            ("Trying to move blank right in rightmost column");
        // move the blank
//        int[][] temp = copy();
//        int blankI = blank/N;
//        int blankJ = blank%N;
//        temp[blankI][blankJ] = (int)tiles[blank+1];
//        temp[blankI][blankJ+1] = 0;
//        Board mv = new Board(temp);
//        return mv;
        int newham = ham;
        int newman = man;
        tiles[blank] = tiles[blank+1];
        tiles[blank+1] = 0;
        if((char) tiles[blank] == blank+2) newham++;
        if((char) tiles[blank] == blank+1) newham--;
        newman -= mandis((blank+1)/N, (blank+1)%N, (char) tiles[blank]);
        newman += mandis(blank/N, blank%N, (char) tiles[blank]);
        Board nb = new Board(tiles, N, newham, newman, blank+1);
        tiles[blank+1] = tiles[blank];
        tiles[blank] = 0;
        return nb;
    }
    
    private static void boardOut(Board x) { // test output
        StdOut.println(x);
        StdOut.print("Dimension: ");
        StdOut.println(x.dimension());
        StdOut.print("Hamming: ");
        StdOut.println(x.hamming());
        StdOut.print("Manhattan: ");
        StdOut.println(x.manhattan());
        StdOut.print("Goal: ");
        StdOut.println(x.isGoal());
        StdOut.println();
    }
        
    public static void main(String[] args){// unit tests (not graded)
        // create initial board from file
        In in = new In(args[0]);
        int N = in.readInt();
        int[][] blocks = new int[N][N];
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++)
            blocks[i][j] = in.readInt();
        Board initial = new Board(blocks);
        
        // print basic functions to standard output
        boardOut(initial);
        
        // print twin
        StdOut.print("Twin: ");
        Board twin = initial.twin();
        boardOut(twin);
                
        // find neighbors and output
        Iterable<Board> nbs = initial.neighbors();
        StdOut.println("Neighbors:");
        for (Board nb : nbs) boardOut(nb);
    }
}