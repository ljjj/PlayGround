# initialization
import numpy
n = 17
a = numpy.zeros(n,dtype = int)
a[:] = -1

# deep first search
def searchThisDigit(a,k,l):
#    if k > 9:
#        print('searching '+str(k)+' of '+str(a))
    if a[n-k] == -1:
        for d in range(10):
            if k==1 and d==0: # first and last digit cannot be 0
                d = 1
            a[n-k] = d
            checkThisDigit(a,k,l)
        a[n-k] = -1 # reset
    else:
        checkThisDigit(a,k,l)
        
def checkThisDigit(a,k,l):
    # find the digit after square using known digits, considering l, the digit added from one order below
    s = a[n-k]*a[n-1]+l
#    if k > 1:
#        s += a[n-k+1]*a[n-2]*10 + a[n-k+1]*a[n-1]
    for i in range(n-k+1,n):
        s += a[i]*a[2*n-k-i-1] # + a[i]*a[2*n-k-i]//10
#    if s % 10 == 9: # possibly missing added 1 from two orders below
#        s = 0
#        for i in range(k):
#            s = s * 10 + a[n - i - 1]
#        print(s)
#        s **= 2
#        for i in range(k - 1):
#            s //= 10
#        s %= 10
#    else:
    l = s // 10 # this adds to the next order digit
    s %= 10
    if a[k-1] == -1:
        a[k-1] = s
        searchThisDigit(a,k+1,l)
        a[k-1] = -1 # reset
    else:
        if a[k-1] == s:
            if k == n:
                print(str(a))
                return
            else:
                checkThisDigit(a,k+1,l)
        else:
            return

searchThisDigit(a,1,0)
