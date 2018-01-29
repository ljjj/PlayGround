/*----------------------------------------------------------------
 *  Author:        Junjiajia Long
 *  Written:       2/13/2015
 *  Last updated:  2/13/2015
 *
 *  Compilation:   javac Subset.java
 *  Execution:     java Subset 3
 *----------------------------------------------------------------*/
public class Subset {
    public static void main(String[] args)
    {
        int k;
        if (args.length < 1) k = 3; // set default k
        else k = Integer.parseInt(args[0]);
        
        RandomizedQueue<String> rq = new RandomizedQueue<String>();
        while (!StdIn.isEmpty()) rq.enqueue(StdIn.readString());
        
        for (int i = 0; i < k; i++) StdOut.println(rq.dequeue());
   }
}