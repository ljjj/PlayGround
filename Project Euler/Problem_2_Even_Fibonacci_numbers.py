from math import log,ceil
from decimal import Decimal
phi = (1+Decimal(5).sqrt())/2
def Fib(n): # return the nth Fibonacci number
    return int(round((phi**n-(-phi)**(-n))/Decimal(5).sqrt()))
def findFib(N): # find the Fibonacci term number not exceeding N
    if N > 21:
        n = int(ceil(Decimal((N*Decimal(5).sqrt())).ln()/Decimal(phi).ln()))
    else:
        n = 8
    while Fib(n) > N:
        n -= 1
    while Fib(n) % 2 != 0:
        n -= 1
    return n
def SumFib(N): # sum over all even Fibonacci number not exceeding N
    n = findFib(N)
    return int(round(((1-phi**(n+1))/(1-phi)-(1-(-1/phi)**(n+1))/(1-(-1/phi)))/Decimal(5).sqrt()))/2
T = input()
for i in range(T):
    print(SumFib(input()))
