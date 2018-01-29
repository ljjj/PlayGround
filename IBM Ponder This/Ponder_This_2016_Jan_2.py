# try to approach this problem with a deep-first search method

# initializations
S = 6
N = 0 # this saves the maximal N found so far and the corresponding sets
NS1 = []
NS2 = []
##paths = [] # saves the path of [S1, S2] searches to avoid repetitions

def matchSets(S1, S2): # given current partial sets S1 and S2 determine the next trial sets
    global N, NS1, NS2 #, paths
##    S1.sort()
##    S2.sort()
##    sets = [S1[:],S2[:]]
##    sets.sort()
##    if sets in paths: # don't repeat
##        return
##    paths.insert(0,sets[:])
    Nm = AllGears(S1[:],S2[:])
    if 0 not in S1 and 0 not in S2: # the sets are complete
        if Nm > N: # replace by the new maximal N and the sets
            N = Nm
            NS1 = [S1[:]]
            NS2 = [S2[:]]
            print(N)
            print(S1)
            print(S2)
            print('')
        elif Nm == N: # register equivalent solutions
            NS1.insert(0,S1[:])
            NS2.insert(0,S2[:])
            print(S1)
            print(S2)
            print('')
    else: # deep-first search
        searchNext(S1[:],S2[:],Nm)
        searchNext(S1[:],S2[:],Nm+1)
        searchNext(S1[:],S2[:],Nm+2)

def searchNext(S1,S2,n): # search factors of n and add to S1 and S2
    cands = factors(n)
    for cand in cands: # cand = [a, b] to be put into the sets
        a = cand[0]
        b = cand[1]
        if a in S1 and 0 in S2:
            newS1 = S1[:]
            newS2 = S2[:]
            newS2[newS2.index(0)] = b
            matchSets(newS1, newS2)
        if b != a and b in S1 and 0 in S2:
            newS1 = S1[:]
            newS2 = S2[:]
            newS2[newS2.index(0)] = a
            matchSets(newS1, newS2)
        if a in S2 and 0 in S1:
            newS1 = S1[:]
            newS2 = S2[:]
            newS1[newS1.index(0)] = b
            matchSets(newS1, newS2)
        if b != a and b in S2 and 0 in S1:
            newS1 = S1[:]
            newS2 = S2[:]
            newS1[newS1.index(0)] = a
            matchSets(newS1, newS2)
        if a not in S1 and b not in S1 and a not in S2 and b not in S2 and 0 in S1 and 0 in S2:
            newS1 = S1[:]
            newS2 = S2[:]
            newS1[newS1.index(0)] = a
            newS2[newS2.index(0)] = b
            matchSets(newS1, newS2)
            if a != b:
                newS1[newS1.index(a)] = b
                newS2[newS2.index(b)] = a
                matchSets(newS1, newS2)

def factors(n): # return a list of factor pair-lists of n
    fs = []
    for i in range(1,int(round(n**0.5))+1):
        if n % i == 0:
            fs.insert(0,[i, n//i])
    return fs

def AllGears(s1,s2): # find the maximal N such that every gear value between 1 and maxN is produced by S1 and S2
    gears = [False] * (max(s1)*max(s2)+2)
    for i in range(0,len(s1)):
        for j in range(0,len(s2)):
            p = s1[i] * s2[j] - 1
            if p > -1:
                gears[p] = True
                gears[p+1] = True
                if p > 0:
                    gears[p-1] = True
    # print(s1,s2,gears)
    return gears.index(False)

S1 = [0]*S
S1[0] = 1
S2 = [0]*S
S2[0] = 1
matchSets(S1[:],S2[:])
S2[0] = 2
matchSets(S1[:],S2[:])
