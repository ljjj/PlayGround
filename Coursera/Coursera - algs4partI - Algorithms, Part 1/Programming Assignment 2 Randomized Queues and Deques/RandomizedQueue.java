/*----------------------------------------------------------------
 *  Author:        Junjiajia Long
 *  Written:       2/13/2015
 *  Last updated:  2/13/2015
 *
 *  Compilation:   javac RandomizedQueue.java
 *  Execution:     java RandomizedQueue
 *----------------------------------------------------------------*/
import java.util.Iterator;

public class RandomizedQueue<Item> implements Iterable<Item> {
    private Item[] rq;
    private int N = 0;
    
    public RandomizedQueue()                 // construct an empty randomized queue
    {
        rq = (Item[]) new Object[1];
    }
    public boolean isEmpty()                 // is the queue empty?
    {
        return N == 0;
    }
    public int size()                        // return the number of items on the queue
    {
        return N;
    }
    public void enqueue(Item item)           // add the item
    {
        if (item == null) throw new java.lang.NullPointerException
            ("Attemp to add a null item");
        if (N == rq.length) resize(2*N);
        rq[N] = item;
        N++;
    }
    public Item dequeue()                    // remove and return a random item
    {
        if (isEmpty()) throw new java.util.NoSuchElementException
            ("Attempt to remove from empty deque");
        int ri = StdRandom.uniform(0,N);
        N--;
        Item item = rq[ri];
        rq[ri] = rq[N];
        rq[N] = null;
        if (N > 0 && N == rq.length/4) resize(2*N);
        return item;
    }
    public Item sample()                     // return (but do not remove) a random item
    {
        if (isEmpty()) throw new java.util.NoSuchElementException
            ("Attempt to remove from empty deque");
        return rq[StdRandom.uniform(0,N)];
    }
    public Iterator<Item> iterator()         // return an independent iterator over items in random order
    {
        return new ListIterator();
    }
    
    private void resize(int capacity)        // resize the item array
    {
        Item[] copy = (Item[]) new Object[capacity];
        for (int i = 0; i < N; i++) copy[i] = rq[i];
        rq = copy;
    }
    private class ListIterator implements Iterator<Item> // implement the iterator
    {
        private Item[] iter = (Item[]) new Object[N];
        private int current = 0;
        
        private ListIterator() // construct the iterator by shuffling
        {
            for (int i = 0; i < N; i++) {
                iter[i] = rq[i];
            }
            for (int i = 0; i < N; i++) {
                int r = StdRandom.uniform(0,i+1);
                Item temp = iter[i];
                iter[i] = iter[r];
                iter[r] = temp;
            }
        }
        
        public boolean hasNext() { return current < N; }
        public void remove() {
            throw new java.lang.UnsupportedOperationException
                ("Attempt to remove in iterator");
        }
        public Item next() {
            if (!hasNext()) throw new java.util.NoSuchElementException
                ("No more items to return in iterator");
            Item item = iter[current];
            iter[current] = null;
            current++;
            return item;
        }
    }
    
    public static void main(String[] args)   // unit testing
    {
        RandomizedQueue<Integer> intRQ = new RandomizedQueue<Integer>();
        StdOut.print("Initial array emptiness: ");
        StdOut.println(intRQ.isEmpty());
        
        StdOut.println("Enqueue 1 3 5 7 9");
        intRQ.enqueue(1);
        intRQ.enqueue(3);
        intRQ.enqueue(5);
        intRQ.enqueue(7);
        intRQ.enqueue(9);
        
        StdOut.print("Emptiness: ");
        StdOut.println(intRQ.isEmpty());
        
        StdOut.print("Sample 10 times: ");
        for (int i = 0; i < 10; i++) StdOut.print(intRQ.sample());
        StdOut.println();
        
        Iterator<Integer> intIter1 = intRQ.iterator();
        StdOut.print("Iterator 1 output: ");
        while (intIter1.hasNext()) StdOut.print(intIter1.next());
        StdOut.println();
        
        Iterator<Integer> intIter2 = intRQ.iterator();
        StdOut.print("Iterator 2 output: ");
        while (intIter2.hasNext()) StdOut.print(intIter2.next());
        StdOut.println();
        
        StdOut.print("Dequeue: ");
        while (!intRQ.isEmpty()) StdOut.print(intRQ.dequeue());
        StdOut.println();
    }
}