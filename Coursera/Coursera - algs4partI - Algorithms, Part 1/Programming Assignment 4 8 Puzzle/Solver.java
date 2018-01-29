/*----------------------------------------------------------------
 *  Author:        Junjiajia Long
 *  Written:       3/7/2015
 *  Last updated:  3/7/2015
 *
 *  Compilation:   javac Solver.java
 *  Execution:     java Solver puzzle00.txt
 *----------------------------------------------------------------*/

public class Solver {
    private final int N; // board dimension < 128
    private final int nSol; // min number of moves to solve initial board
    private final Stack<Board> sol; // solution
    
    public Solver(Board initial) {         // find a solution to the initial board (using the A* algorithm)
        if (initial == null)
            throw new NullPointerException("Initial board null");
        N = initial.dimension();        
//        // initialize the game and its twin
//        Node currO = new Node(initial,null); // root node
//        Node prevO = currO; // previous search node
//        MinPQ<Node> gameTreeO = new MinPQ<Node>(); // game tree of the original board
//        Node currT = new Node(initial.twin(),null); // root node
//        Node prevT = currT; // previous search node
//        MinPQ<Node> gameTreeT = new MinPQ<Node>(); // game tree of the twin board
//        
//        // A* algorithm simultaneously on both original and twin board
//        while (!currO.board.isGoal() && !currT.board.isGoal()) { // run until goal board is found in either tree
//            Iterable<Board> nbsO = currO.board.neighbors(); // neighbors of current node
//            for (Board nb : nbsO)
//                if (nb != prevO.board)
//                    gameTreeO.insert(new Node(nb,currO));
//            prevO = currO;
//            currO = gameTreeO.delMin(); // least priority node
//            
//            Iterable<Board> nbsT = currT.board.neighbors(); // neighbors of current node
//            for (Board nb : nbsT)
//                if (nb != prevT.board) // critical optimization
//                    gameTreeT.insert(new Node(nb,currT));
//            prevT = currT;
//            currT = gameTreeT.delMin(); // least priority node
//        }
//        
//        // define the solution
//        sol = new Stack<Board>();
//        if (currO.board.isGoal()) {
//            nSol = currO.nMoves;
//            while (currO != null) {
//                sol.push(currO.board);
//                currO =  currO.prev;
//            }
//        } else {
//            nSol = -1;
//        }
        sol = new Stack<Board>();
        checkParity cP = new checkParity(initial);
        if (!cP.solvable) nSol = -1;
        else {
            // initialize the game and its twin
            Node curr = new Node(initial,null,cP.blank); // root node
            MinPQ<Node> gameTree = new MinPQ<Node>(); // game tree
            // A* algorithm
            while (!curr.board.isGoal()) { // run until goal board is found
                Iterable<Board> nbs = curr.board.neighbors(); // neighbors of current node
                // from current blank determine neighbors' blank
                int[] newBlank;
                if (curr.blank == 0) { // top-left corner
                    newBlank = new int[2];
                    newBlank[0] = 1;
                    newBlank[1] = N;
                } else if (curr.blank == N-1) { // top-right corner
                    newBlank = new int[2];
                    newBlank[0] = N-2;
                    newBlank[1] = 2*N-1;
                } else if (curr.blank == N*(N-1)) { // bottom-left corner
                    newBlank = new int[2];
                    newBlank[0] = N*(N-1)+1;
                    newBlank[1] = N*(N-2);
                } else if (curr.blank == N*N-1) { // bottom-right corner
                    newBlank = new int[2];
                    newBlank[0] = N*N-2;
                    newBlank[1] = N*(N-1)-1;
                } else if (curr.blank < N) { // top row non-corner
                    newBlank = new int[3];
                    newBlank[0] = curr.blank+1;
                    newBlank[1] = curr.blank-1;
                    newBlank[2] = curr.blank+N;
                } else if (curr.blank > N*(N-1)) { // bottom row non-corner
                    newBlank = new int[3];
                    newBlank[0] = curr.blank+1;
                    newBlank[1] = curr.blank-1;
                    newBlank[2] = curr.blank-N;
                } else if (curr.blank % N == 0) { // leftmost column non-corner
                    newBlank = new int[3];
                    newBlank[0] = curr.blank+1;
                    newBlank[1] = curr.blank+N;
                    newBlank[2] = curr.blank-N;
                } else if (curr.blank % N == N-1) { // rightmost column non-corner
                    newBlank = new int[3];
                    newBlank[0] = curr.blank-1;
                    newBlank[1] = curr.blank+N;
                    newBlank[2] = curr.blank-N;
                } else { // non-edge
                    newBlank = new int[4];
                    newBlank[0] = curr.blank+1;
                    newBlank[1] = curr.blank-1;
                    newBlank[2] = curr.blank+N;
                    newBlank[3] = curr.blank-N;
                }
                
                // check repetition by blank equality
                int nB = 0;
                for (Board nb : nbs) {
//                    boolean foundPrev = false;
//                    if ( !foundPrev && (curr.nMoves == 0 || nb.hamming() != curr.prev.board.hamming()
//                    || nb.manhattan() != curr.prev.board.manhattan()
//                    || nb.m != curr.prev.board) ) {
                    if (curr.nMoves == 0)
                        gameTree.insert(new Node(nb,curr,newBlank[nB]));
                    else if (newBlank[nB] != curr.prev.blank)
                        gameTree.insert(new Node(nb,curr,newBlank[nB]));
//                    } else foundPrev = true;
                    nB++;
                }
//                StdOut.print(curr.priority);
//                StdOut.print(" ");
//                StdOut.println(curr.nMoves);
                curr = gameTree.delMin(); // least priority node
            }
            
            // get the solution
            nSol = curr.nMoves;
            while (curr != null) {
                sol.push(curr.board);
                curr =  curr.prev;
            }
        }
    }
    
    private final class checkParity { // check the parity and blank position of the board
        public final boolean solvable; // whether the board is solvable
        public final int blank; // position of the blank
        
        public checkParity(Board b) {
            // initalize
            int parity = 0;
            String[] s = b.toString().split("\n"); // split the board string into rows
            char[] temp = new char[N*N];
            int bla = -1;
            
            // split the row string into individual keys
            for (int i = 0; i < N; i++) {
                String r[] = s[i+1].split(" ");
                int c = 0;
                for (int j = 0; j < N; j++) {
                    while (r[c].equals("")) c++;
                    temp[N*i+j] = (char) Integer.parseInt(r[c]);
                    if (temp[N*i+j] == 0) bla = N*i+j;
                    c++;
                }
            }
            blank = bla;
            assert blank >= 0;
            for (int k = 0; k < N*N; k++)
                if (k != blank)
                for (int l = k+1; l < N*N; l++)
                if (l != blank)
                if (temp[k] > temp[l]) parity++;
            if (N % 2 == 0) // even-sized board need to add the row number of blank
                parity = parity + blank/N + 1;
            if (parity % 2 == 0) solvable = true;
            else solvable = false;
        }
    }
    
    public boolean isSolvable() {          // is the initial board solvable?
        return nSol != -1;
    }
    public int moves() {                   // min number of moves to solve initial board; -1 if unsolvable
        return nSol;
    }
    public Iterable<Board> solution() {    // sequence of boards in a shortest solution; null if unsolvable
        if (!isSolvable()) return null;
        return sol;
    }
    
    private class Node implements Comparable<Node>{ // search node
        public final Board board; // configuration of the node
        public final int nMoves; // number of moves made to reach the node
        public final int priority; // priority of this node
        public final Node prev; // previous search node
        public final int blank; // blank position
        
        public Node(Board c, Node p, int b) { // construct the search node based on the previous node
            board = c;
            if (p != null) {
                nMoves = p.nMoves + 1;
                prev = p;}
            else {
                nMoves = 0;
                prev = null;
            }
            blank = b;
            priority = nMoves + board.manhattan();
        }
        
        public int compareTo(Node that) { // does this node have smaller priority?
            if (this.priority < that.priority) { return -1;}
            else if (this.priority > that.priority) { return 1;}
//            else if (this.board.hamming() < that.board.hamming()) { return -1;}
//            else if (this.board.hamming() > that.board.hamming()) { return 1;}
            else if (this.board.manhattan() < that.board.manhattan()) { return -1;}
            else if (this.board.manhattan() > that.board.manhattan()) { return 1;}
            { return 0;}
        }
    }
    
    public static void main(String[] args){// solve a slider puzzle (given below)
        // create initial board from file
        In in = new In(args[0]);
        int N = in.readInt();
        int[][] blocks = new int[N][N];
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++)
            blocks[i][j] = in.readInt();
        Board initial = new Board(blocks);
        
        // solve the puzzle
        Solver solver = new Solver(initial);
        
        // print solution to standard output
        if (!solver.isSolvable())
            StdOut.println("No solution possible");
        else {
            StdOut.println("Minimum number of moves = " + solver.moves());
            for (Board board : solver.solution())
                StdOut.println(board);
        }
    }
}