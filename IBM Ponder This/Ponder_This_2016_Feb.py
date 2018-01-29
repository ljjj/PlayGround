# Step 1: we look at n**i = ban, the only possibilities are
# 2**9=512 eliminated by Step 2
# 3**5=243
# 5**3=125
# 5**4=625
# 6**3=216 eliminated by Step 2
# 9**3=729
# Step 2: since 1**x = 1, we eliminate the possibility that a = 1.
# Step 3: c != 0 and t != 0 as they are the beginning of a number
# Step 3: brute force the remaining possibilities
from itertools import permutations

def checkEqs(a,b,c,e,i,m,n,o,r,t):
    connect = (((((c*10+o)*10+n)*10+n)*10+e)*10+c)*10+t
    inter = (((i*10+n)*10+t)*10+e)*10+r
    toe = (t*10+o)*10+e
    at = a*10+t
    cm = c*10+m
    aconnect = a**connect
    atoe = a**toe
    return aconnect % inter == c % inter and aconnect % at == atoe % at and aconnect % cm == atoe % cm

def checkRemains(n,i,b,a,remains):
    for [c,e,m,o,r,t] in permutations(remains):
        if c != 0 and t != 0:
            if checkEqs(a,b,c,e,i,m,n,o,r,t):
                print('a:'+str(a))
                print('b:'+str(b))
                print('c:'+str(c))
                print('e:'+str(e))
                print('i:'+str(i))
                print('m:'+str(m))
                print('n:'+str(n))
                print('o:'+str(o))
                print('r:'+str(r))
                print('t:'+str(t))

n,i,b,a = 3,5,2,4
remains = [0,1,6,7,8,9]
checkRemains(n,i,b,a,remains)
n,i,b,a = 5,3,1,2
remains = [0,4,6,7,8,9]
checkRemains(n,i,b,a,remains)
n,i,b,a = 5,4,6,2
remains = [0,1,3,7,8,9]
checkRemains(n,i,b,a,remains)
n,i,b,a = 9,3,7,2
remains = [0,1,4,5,6,8]
checkRemains(n,i,b,a,remains)
