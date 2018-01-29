/*----------------------------------------------------------------
 *  Author:        Junjiajia Long
 *  Written:       1/25/2015
 *  Last updated:  1/30/2015
 *
 *  Compilation:   javac PercolationStats.java
 *  Execution:     java PercolationStats 200 100
 *----------------------------------------------------------------*/
public class PercolationStats {
    private int T;
    private double m; // calculate mean
    private double std; // calculate stddev
    
    public PercolationStats(int N, int T)
    {// perform T independent experiments on an N-by-N grid
        if (N <= 0)
            throw new IllegalArgumentException("Grid size N must be positive");
        if (T <= 0)
            throw new IllegalArgumentException
            ("Number of independent computations T must be positive");
        this.T = T;
        Percolation[] experiments = new Percolation[T];
        double[] x = new double[T];
        for (int t = 0; t < T; t++) // perform each experiment
        {
            experiments[t] = new Percolation(N);
            x[t] = 0;
            while (!experiments[t].percolates()) // open site until percolation
            {
                // choose a site randomly and open
                int i = StdRandom.uniform(1,N+1);
                int j = StdRandom.uniform(1,N+1);
                //StdOut.print(i);
                //StdOut.print(" ");
                //StdOut.println(j);
                // be careful not to double count
                if (!experiments[t].isOpen(i,j))
                {
                    experiments[t].open(i,j);
                    x[t]++; // add to total number of open sites
                }
            }
            x[t] /= (N*N);
        }
        m = StdStats.mean(x);
        std = StdStats.stddev(x);
    }
    
    public double mean() // sample mean of percolation threshold
    {
        return m;
    }
    public double stddev() // sample standard deviation of percolation threshold
    {
        return std;
    }
    public double confidenceLo() // low  endpoint of 95% confidence interval
    {
        return m - 1.96*std/Math.sqrt(T);
    }
    public double confidenceHi() // high endpoint of 95% confidence interval
    {
        return m + 1.96*std/Math.sqrt(T);
    }
    
    public static void main(String[] args)    // test client (described below)
    {
        int T;
        int N;
        if (args.length < 2) T = 100; // set default T
        else T = Integer.parseInt(args[1]);
        if (args.length < 1) N = 200; // set default N
        else N = Integer.parseInt(args[0]);
        // perform computational experiments
        PercolationStats pStats = new PercolationStats(N, T);
        // output stats
        StdOut.print("mean                    = ");
        StdOut.println(pStats.mean());
        StdOut.print("stddev                  = ");
        StdOut.println(pStats.stddev());
        StdOut.print("95% confidence interval = ");
        StdOut.print(pStats.confidenceLo());
        StdOut.print(", ");
        StdOut.println(pStats.confidenceHi());
    }
}