from math import sqrt
def factors(n): # return number of factors
    s = 0
    for i in range(int(sqrt(n))):
        if n%(i+1) == 0:
            s += 2
    if sqrt(n)*sqrt(n) == n:
        s -= 1
    return s

k = 1
while factors(k*(k+1)/2) <= 500:
    k += 1
print(k*(k+1)/2)
print(factors(k*(k+1)/2))
