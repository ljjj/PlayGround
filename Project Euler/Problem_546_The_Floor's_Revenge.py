from math import floor
import numpy

N = 1
K = 50
kOffset = 4
## method 1: fs = numpy.zeros((N,K),dtype=long)

def f(k,n):
##    global fs
    if n < k:
        return n+1
## method 1:    if fs[k-kOffset][n] != 0:
## method 1:        return fs[k-kOffset][n]
## method 1:    s = 0
## method 1:    for i in range(int(n)+1):
## method 1:        s += f(k,floor(i/k))
## method 1:    fs[k-kOffset][n] = s
## method 1:    return s
    else:
## method 2:        return f(k,n-n%k-1)+(n%k+1)*f(k,n/k)
        s = 0
        for i in range(n/k):
            s += k*f(k,i) % 1000000007
        return s + (n%k+1)*f(k,n/k) % 1000000007

    
for i in range(N):
    for j in range(K):
        print('f_'+str(i+kOffset)+'('+str(j)+')='+str(f(i+kOffset,j)))
