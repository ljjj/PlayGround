os = 2
Nm = 1000000 - os
primes = []
sums = []
for i in range(Nm+1):
    primes.append(True)
    sums.append(0)
s = 0
for i in range(Nm+1):
    if primes[i]:
        s += i+os
        sums[i] = s
        x = i+os
        y = x+x
        while y-os <= Nm:
            primes[y-os] = False
            y += x

def sumPrimes(N):
    j = N-os
    while not primes[j]:
        j = j-1
    return sums[j]

T = input()
for i in range(T):
    print(sumPrimes(input()))
