import numpy as np
import pickle as pk

# define possible movements, p is a 2-tuple
def Qmove(p):
    move = np.zeros((n,n), dtype=bool)
    for i in range(n):
        for j in range(n):
            if i == p[0] or j == p[1] or abs(p[0]-i) == abs(p[1]-j):
                move[i,j] = True
    return move

def Kmove(p):
    move = np.zeros((n,n), dtype=bool)
    for i in range(n):
        for j in range(n):
            if abs(p[0]-i) <= 1 and abs(p[1]-j) <= 1:
                move[i,j] = True
    return move

def Rmove(p):
    move = np.zeros((n,n), dtype=bool)
    for i in range(n):
        for j in range(n):
            if i == p[0] or j == p[1]:
                move[i,j] = True
    return move

def Bmove(p):
    move = np.zeros((n,n), dtype=bool)
    for i in range(n):
        for j in range(n):
            if abs(p[0]-i) == abs(p[1]-j):
                move[i,j] = True
    return move

def Nmove(p):
    move = np.zeros((n,n), dtype=bool)
    for i in range(n):
        for j in range(n):
            if abs((p[0]-i) * (p[1]-j)) == 2:
                move[i,j] = True
    return move


# find the moves and loop over depth first search
def findMove(m, chessboard, Q, K, R, B, N):
    # determine avaiable moves
    Qnext = np.where(np.logical_and(chessboard == -1, Qmove(Q[-1])))
    Knext = np.where(np.logical_and(chessboard == -1, Kmove(K[-1])))
    Rnext = np.where(np.logical_and(chessboard == -1, Rmove(R[-1])))
    Bnext = np.where(np.logical_and(chessboard == -1, Bmove(B[-1])))
    Nnext = np.where(np.logical_and(chessboard == -1, Nmove(N[-1])))
    
    # loop into possible moves
    for qi in range(len(Qnext[0])):
        qnext = (Qnext[0][qi], Qnext[1][qi])
        # find new products for this move
        chessboard[qnext] = m # tentative move
        rp = chessboard[qnext[0],:]
        rp = np.prod(rp[np.where(rp > 0)] % 10)
        cp = chessboard[:,qnext[1]]
        cp = np.prod(cp[np.where(cp > 0)] % 10)
        if rp == 0 or cp == 0:
            print(chessboard)
            print(qnext)
            print(m)
        # check new products are factors
        if rprod[qnext[0]] % rp != 0 or cprod[qnext[1]] % cp != 0:
            chessboard[qnext] = -1 # remove tentative move
            continue
        # places under attack
        attacked_q = Qmove(qnext)
        
        for ki in range(len(Knext[0])):
            knext = (Knext[0][ki], Knext[1][ki])
            # here is under attack
            if attacked_q[knext]:
                continue
            # here attack others
            attacked_k = Kmove(knext)
            if attacked_k[qnext]:
                continue
            # find new products for this move
            chessboard[knext] = m + 10 # tentative move
            rp = chessboard[knext[0],:]
            rp = np.prod(rp[np.where(rp > 0)] % 10)
            cp = chessboard[:,knext[1]]
            cp = np.prod(cp[np.where(cp > 0)] % 10)
            # check new products are factors
            if rprod[knext[0]] % rp != 0 or cprod[knext[1]] % cp != 0:
                chessboard[knext] = -1 # remove tentative move
                continue
            
            for ri in range(len(Rnext[0])):
                rnext = (Rnext[0][ri], Rnext[1][ri])
                # here is under attack
                if attacked_q[rnext] or attacked_k[rnext]:
                    continue
                # here attack others
                attacked_r = Rmove(rnext)
                if attacked_r[qnext] or attacked_r[knext]:
                    continue
                # find new products for this move
                chessboard[rnext] = m + 20 # tentative move
                rp = chessboard[rnext[0],:]
                rp = np.prod(rp[np.where(rp > 0)] % 10)
                cp = chessboard[:,rnext[1]]
                cp = np.prod(cp[np.where(cp > 0)] % 10)
                # check new products are factors
                if rprod[rnext[0]] % rp != 0 or cprod[rnext[1]] % cp != 0:
                    chessboard[rnext] = -1 # remove tentative move
                    continue
                
                for bi in range(len(Bnext[0])):
                    bnext = (Bnext[0][bi], Bnext[1][bi])
                    # here is under attack
                    if attacked_q[bnext] or attacked_k[bnext] or attacked_r[bnext]:
                        continue
                    # here attack others
                    attacked_b = Bmove(bnext)
                    if attacked_b[qnext] or attacked_b[knext] or attacked_b[rnext]:
                        continue
                    # find new products for this move
                    chessboard[bnext] = m + 30 # tentative move
                    rp = chessboard[bnext[0],:]
                    rp = np.prod(rp[np.where(rp > 0)] % 10)
                    cp = chessboard[:,bnext[1]]
                    cp = np.prod(cp[np.where(cp > 0)] % 10)
                    # check new products are factors
                    if rprod[bnext[0]] % rp != 0 or cprod[bnext[1]] % cp != 0:
                        chessboard[bnext] = -1 # remove tentative move
                        continue
                    
                    for ni in range(len(Nnext[0])):
                        nnext = (Nnext[0][ni], Nnext[1][ni])
                        # here is under attack
                        if attacked_q[nnext] or attacked_k[nnext] or attacked_r[nnext] or attacked_b[nnext]:
                            continue
                        # here attack others
                        attacked_n = Nmove(nnext)
                        if attacked_n[qnext] or attacked_n[knext] or attacked_n[rnext] or attacked_n[bnext]:
                            continue
                        # find new products for this move
                        chessboard[nnext] = m + 40 # tentative move
                        rp = chessboard[nnext[0],:]
                        rp = np.prod(rp[np.where(rp > 0)] % 10)
                        cp = chessboard[:,nnext[1]]
                        cp = np.prod(cp[np.where(cp > 0)] % 10)
                        # check new products are factors
                        if rprod[nnext[0]] % rp != 0 or cprod[nnext[1]] % cp != 0:
                            chessboard[nnext] = -1 # remove tentative move
                            continue
                        
                        moveIt(m, chessboard, Q, K, R, B, N,
                               qnext, knext, rnext, bnext, nnext)
                        
                        chessboard[nnext] = -1
                    chessboard[bnext] = -1
                chessboard[rnext] = -1
            chessboard[knext] = -1
        chessboard[qnext] = -1


# confirm moves
def moveIt(m, chessboard, Q, K, R, B, N,
           qnext, knext, rnext, bnext, nnext):
    global largestpsum, optimalboard
    
    # push to stack
    Q.append(qnext)
    K.append(knext)
    R.append(rnext)
    B.append(bnext)
    N.append(nnext)
    
    # find products
    rprodcurr, cprodcurr = boardProd(chessboard, Q, K, R, B, N)
    
    if m > 5:
        print(rprodcurr, cprodcurr)
        print(chessboard)
    
    # start searching the next depth
    if m < 8:
        if m < 7: # check product before 7 moves as factors
            checkProd = (np.remainder(rprod, rprodcurr) == 0).all() and (np.remainder(cprod, cprodcurr) == 0).all()
        elif m == 7: # check product after 7 moves to be equal
            checkProd = (rprodcurr == rprod).all() and (cprodcurr == cprod).all()
        if checkProd:
            findMove(m + 1, chessboard, Q, K, R, B, N)
    
    # check final product sum
    elif m == 8:
        print(chessboard)
        psum = sum(rprodcurr) + sum(cprodcurr)
        if psum > largestpsum:
            largestpsum = psum
            optimalboard = chessboard.copy()
            print(largestpsum)
            print(chessboard)
            print('')
    
    # pop from stack
    Q.pop()
    K.pop()
    R.pop()
    B.pop()
    N.pop()


# calculate board products
def boardProd(chessboard, Q, K, R, B, N):
    # remove 0
    chessboard[Q[0]] = 1
    chessboard[K[0]] = 1
    chessboard[R[0]] = 1
    chessboard[B[0]] = 1
    chessboard[N[0]] = 1
    
    # remove -1
    negind = np.where(chessboard < 0)
    chessboard[negind] = 1
    
    # -1 is taken care of by abs
    rprodcurr = abs(np.prod(chessboard % 10, axis=1))
    cprodcurr = abs(np.prod(chessboard % 10, axis=0))
    
    # put back 0 and -1
    chessboard[Q[0]] = 0
    chessboard[K[0]] = 0
    chessboard[R[0]] = 0
    chessboard[B[0]] = 0
    chessboard[N[0]] = 0
    chessboard[negind] = -1
    
    return rprodcurr, cprodcurr


# initialization
n = 8 # board size
rprod = np.array([7, 1890, 8, 10080, 20, 840, 144, 1260]) # row products after 7 moves
cprod = np.array([2744, 36, 375, 336, 108, 240, 20, 504]) # column products after 7 moves
chessboard = np.ones((n,n), dtype=int)*-1 # board
largestpsum = 0
optimalboard = []
# stack of positions
Q = []
Q.append((0,1))
K = []
K.append((1,7))
R = []
R.append((4,4))
B = []
B.append((7,6))
N = []
N.append((6,0))
# first push into stack
chessboard[Q[-1]] = 0
chessboard[K[-1]] = 0
chessboard[R[-1]] = 0
chessboard[B[-1]] = 0
chessboard[N[-1]] = 0

# depth first search
findMove(1, chessboard, Q, K, R, B, N)

# save results
with open('Puzzle_2016_Sep_result.pckl','wb') as f:
    pk.dump(largestpsum, f)
    pk.dump(optimalboard, f)
with open('Puzzle_2016_Sep_result.txt','w') as f:
    f.write(str(largestpsum))
    f.write('\n')
    f.write(str(optimalboard))