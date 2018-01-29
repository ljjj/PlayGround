/*----------------------------------------------------------------
 *  Author:        Junjiajia Long
 *  Written:       1/25/2015
 *  Last updated:  1/29/2015
 *
 *  Compilation:   javac Percolation.java
 *  Execution:     java Percolation
 *----------------------------------------------------------------*/
public class Percolation {
    private int N;
    // the acutal sites are defined in a 1D array with index (i-1)*N+j
    // the top virtual site is identified with 0 and the bottom N^2+1
    private boolean[] sites; // open or blocked state of the site
    private WeightedQuickUnionUF fullUF; // used to determine if full
    private WeightedQuickUnionUF percUF; // used to determine if percolates
    
    public Percolation(int N) // create N-by-N grid, with all sites blocked
    {
        if (N <= 0)
            throw new IllegalArgumentException("Grid size N must be positive");
        this.N = N;
        sites = new boolean[N*N+2];
        fullUF = new WeightedQuickUnionUF(N*N+2);
        percUF = new WeightedQuickUnionUF(N*N+2);
        // initialize virtual sites
        sites[0] = true;
        sites[N*N+1] = true;
        // initialize actual sites
        for (int i = 1; i <= N; i++)
        {
            for (int j = 1; j <= N; j++)
            {
                sites[(i-1)*N+j] = false;
            }
        }
    }

    public void open(int i, int j) // open site (row i, column j)
    { // if it is not open already
        checkIndex(i,j);
        if (!isOpen(i,j))
        {
            sites[(i-1)*N+j] = true;
            if (i == 1) // connect to virtual top
            {
                fullUF.union((i-1)*N+j, 0);
                percUF.union((i-1)*N+j, 0);
            }
            else if (sites[(i-2)*N+j]) // connect to site above if it is open
            {
                fullUF.union((i-1)*N+j, (i-2)*N+j);
                percUF.union((i-1)*N+j, (i-2)*N+j);
            }
            if (i == N) // only percUF connects to virtual bottom
            {
                percUF.union((i-1)*N+j, N*N+1);
            }
            else if (sites[i*N+j]) // connect to site below if it is open
            {
                fullUF.union((i-1)*N+j, i*N+j);
                percUF.union((i-1)*N+j, i*N+j);
            }
            if ((j > 1) && (sites[(i-1)*N+j-1])) // connect to left open site
            {
                fullUF.union((i-1)*N+j, (i-1)*N+j-1); 
                percUF.union((i-1)*N+j, (i-1)*N+j-1);
            }
            if ((j < N) && (sites[(i-1)*N+j+1])) // connect to right open site
            {
                fullUF.union((i-1)*N+j, (i-1)*N+j+1);
                percUF.union((i-1)*N+j, (i-1)*N+j+1); 
            }
        }
    }
    public boolean isOpen(int i, int j) // is site (row i, column j) open?
    {
        checkIndex(i,j);
        return sites[(i-1)*N+j];
    }
    public boolean isFull(int i, int j) // is site (row i, column j) full?
    {
        checkIndex(i,j);
        return fullUF.connected((i-1)*N+j, 0);
    }
    public boolean percolates() // does the system percolate?
    {
        return percUF.connected(0, N*N+1);
    }
    
    private void checkIndex(int i, int j) // check if index is out of bounds
    {
        if (i <= 0 || i > N)
            throw new IndexOutOfBoundsException("Row index i out of bounds");
        if (j <= 0 || j > N)
            throw new IndexOutOfBoundsException("Column index j out of bounds");
    }
    private void testLattice() // output the site values and percolation state
    { // used for testing
        for (int i = 1; i <= N; i++)
        {
            for (int j = 1; j <= N; j++)
                if (isFull(i,j)) StdOut.print(2);
                else if (isOpen(i,j)) StdOut.print(1);
                else StdOut.print(0);
            StdOut.println();
        }
        StdOut.println(percolates());
    }

    public static void main(String[] args){// test client (optional)
        int N = 4;
        Percolation perc = new Percolation(N);
        // construct a lattice and test it
        perc.open(1,1);
        perc.open(2,1);
        perc.open(2,4);
        perc.open(4,1);
        perc.open(3,2);
        perc.open(3,4);
        perc.testLattice();
        perc.open(3,1);
        perc.open(3,3);
        perc.testLattice();
    }
}