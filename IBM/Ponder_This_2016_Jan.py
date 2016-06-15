import sys, time

# initializations
S = 6
N = 0 # this saves the maximal N found so far and the corresponding sets
NS1 = []
NS2 = []

def testMaxWithSetSize(M, size1, size2): # find maximal N with size1 size2 sets constrained to be <= M
    global N
    global NS1
    global NS2
    if M < size1 or M < size2: # the numbers in the sets should all be distinct, same for S2
        sys.stderr.write('M < size1 or M < size2 in testMaxSetWithSetSize')
        sys.exit(1)
    si1 = list(range(M-size1+1,M+1))
    si1.reverse()
    # loop over all possible sets of S1 keeping S1[0] = M
    while si1[0] == M:
        # loop over all possible sets of S2 starting with S1
        if size2 == size1: # exploit symmetry
            si2 = si1[:]
        else:
            si2 = list(range(M-size2+1,M+1))
            si2.reverse()
        while True:
            Ncurr = AllGears(si1,si2)
            if Ncurr > N: # replace by the new maximal N and the sets
                N = Ncurr
                NS1 = [si1[:]]
                NS2 = [si2[:]]
                print(N)
                print(NS1)
                print(NS2)
            elif Ncurr == N: # register equivalent solutions
                NS1.insert(0,si1)
                NS2.insert(0,si2)
            if si2[0] == size2:
                break
            si2 = advanceSet(si2)
        if si1[0] == size1:
            break
        si1 = advanceSet(si1)

def advanceSet(Si): # find the next set
    Ns = len(Si)
    if Si[0] == Ns: # assuming Si is sorted in reverse order
        sys.stderr.write('cannot advance')
        sys.exit(1)
    k = Ns-1
    l = 1
    while Si[k] == l: # find the number to change
        k = k - 1 # k >= 0 or else the above error catches it
        l = l + 1
    Si[k] = Si[k] - 1 # all numbers after it are consecutive decreasing integers
    temp = list(range(Si[k]-Ns+k+1,Si[k]))
    temp.reverse()
    Si[(k+1):Ns] = temp[:]
    return Si

def AllGears(s1,s2): # find the maximal N such that every gear value between 1 and maxN is produced by S1 and S2
    gears = [False] * (max(s1)*max(s2)+2)
    for i in range(0,len(s1)):
        for j in range(0,len(s2)):
            p = s1[i] * s2[j] - 1
            gears[p] = True
            gears[p+1] = True
            if p > 0:
                gears[p-1] = True
    # print(s1,s2,gears)
    return gears.index(False)


# there are a total of S^2 products and each covers at most 3 gears, so N
# will not exceed 3*S^2. If the maximum of S1 and S2 exceeds that it is useless
for m in range(S,3*S**2+1):
    print('Testing m = '+str(m))
    start_time = time.time()
    testMaxWithSetSize(m,S,S)
    elapsed_time = time.time() - start_time
    print('Elapsed time is '+str(elapsed_time))
