from math import sqrt,log,ceil
# precalculate relevant Fibonacci numbers
phi = (1+sqrt(5))/2
k = int(ceil(log(40000000000000000*sqrt(5))/log(phi)))
Fibs = []
for i in range(k):
    Fibs.append(0)
Fibs[0] = 1
Fibs[1] = 2
def Fib(n): # return the nth Fibonacci number
    global Fibs
    if Fibs[n] == 0:
        Fibs[n] = Fib(n-1)+Fib(n-2)
    return Fibs[n]
def findFib(N): # find the Fibonacci term number not exceeding N
    n = 0
    while Fib(n) <= N:
       n += 1
    return n-1
def SumFib(N): # sum over all even Fibonacci number not exceeding N
    n = findFib(N)
    s = 0
    for i in range(n):
       if Fibs[i] % 2 == 0:
           s += Fibs[i]
    return s
T = input()
for i in range(T):
    print(SumFib(input()))
