import numpy
os = 2
N = 2000000-os
prime = numpy.ones(N, dtype = bool)
for i in range(N):
    if prime[i]:
        x = i+os
        y = x+x
        while y-os < N:
            prime[y-os] = False
            y += x
s = 0
for i in range(N):
    if prime[i]:
        x = i + os
        s += x
print(s)
