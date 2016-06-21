# determine delegate assignments
def delegates(xs, n, d):
    cs = [d*x//n for x in xs]
    os = [float(d)*x/n-d*x//n for x in xs]
    os = sorted(range(len(os)), key=lambda k: os[k])
    cremain = d - sum(cs)
    j = len(os)
    while cremain > 0:
        j -= 1
        cs[os[j]] += 1
        cremain -= 1
    return cs

# initial counting
N = 102
D = 7
# loop over possible ballots cast
for xH in range(N+1):
    for xL in range(N-xH+1):
        xM = N - xH - xL
        H, L, M = delegates([xH,xL,xM], N, D)
        if H == 2 and L == 2 and M == 3:
            # after audit
            H, L, M = delegates([xH+1,xL,xM], N+1, D+1)
            if H == 1:
                print xH, xL, xM