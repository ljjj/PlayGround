/*----------------------------------------------------------------
 *  Author:        Junjiajia Long
 *  Written:       2/12/2015
 *  Last updated:  2/13/2015
 *
 *  Compilation:   javac Deque.java
 *  Execution:     java Deque
 *----------------------------------------------------------------*/
import java.util.Iterator;

public class Deque<Item> implements Iterable<Item> {
    private Node first = null;
    private Node last = null;
    private int N = 0; // deque size
    private class Node
    {
        Item item;
        Node prev;
        Node next;
    }

    public Deque()                           // construct an empty deque
    {
    }
    public boolean isEmpty()                 // is the deque empty?
    {
        return N == 0;
    }
    public int size()                        // return the number of items on the deque
    {
        return N;
    }
    public void addFirst(Item item)          // add the item to the front
    {
        if (item == null) throw new java.lang.NullPointerException
            ("Attemp to add a null item");
        Node oldfirst = first;
        first = new Node();
        first.item = item;
        first.prev = null;
        first.next = oldfirst;
        if (isEmpty()) last = first;
        else  oldfirst.prev = first;
        N++;
    }
    public void addLast(Item item)           // add the item to the end
    {
        if (item == null) throw new java.lang.NullPointerException
            ("Attemp to add a null item");
        Node oldlast = last;
        last = new Node();
        last.item = item;
        last.prev = oldlast;
        last.next = null;
        if (isEmpty()) first = last;
        else oldlast.next = last;
        N++;
    }
    public Item removeFirst()                // remove and return the item from the front
    {
        if (isEmpty()) throw new java.util.NoSuchElementException
            ("Attempt to remove from empty deque");
        Item item = first.item;
        first = first.next;
        N--;
        if (isEmpty()) last = null;
        else     first.prev = null;
        return item;
    }
    public Item removeLast()                 // remove and return the item from the end
    {
        if (isEmpty()) throw new java.util.NoSuchElementException
            ("Attempt to remove from empty deque");
        Item item = last.item;
        last = last.prev;
        N--;
        if (isEmpty()) first = null;
        else       last.next = null;
        return item;
    }
    public Iterator<Item> iterator()         // return an iterator over items in order from front to end
    {
        return new ListIterator();
    }
   
    private class ListIterator implements Iterator<Item> // implement the iterator
    {
        private Node current = first;
        
        public boolean hasNext() { return current != null; }
        public void remove() {
            throw new java.lang.UnsupportedOperationException
                ("Attempt to remove in iterator");
        }
        public Item next() {
            if (!hasNext()) throw new java.util.NoSuchElementException
                ("No more items to return in iterator");
            Item item = current.item;
            current = current.next;
            return item;
        }
    }
    
    public static void main(String[] args)   // unit testing
    {
        Deque<Integer> intDeq = new Deque<Integer>();
        StdOut.print("Initial array emptiness: ");
        StdOut.println(intDeq.isEmpty());
        
        StdOut.println("Add 1 3 5 first, add 7 9 last");
        intDeq.addFirst(1);
        intDeq.addFirst(3);
        intDeq.addFirst(5);
        intDeq.addLast(7);
        intDeq.addLast(9);
        
        StdOut.print("Emptiness: ");
        StdOut.println(intDeq.isEmpty());
        
        Iterator<Integer> intIter = intDeq.iterator();
        StdOut.print("Iterator output: ");
        while (intIter.hasNext()) StdOut.print(intIter.next());
        StdOut.println();
        
        StdOut.println("Remove first once, then remove last till empty");
        StdOut.print("Pop order: ");
        StdOut.print(intDeq.removeFirst());
        while (!intDeq.isEmpty()) StdOut.print(intDeq.removeLast());
        StdOut.println();
    }
}