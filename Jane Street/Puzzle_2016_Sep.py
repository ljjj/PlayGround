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
            if abs(p[0]-i) * abs(p[1]-j) == 2:
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
        attacked = Qmove(qnext) # places under attack
        
        for ki in range(len(Knext[0])):
            knext = (Knext[0][ki], Knext[1][ki])
            if attacked[knext]: # cannot move here
                continue
            attacked = np.logical_and(attacked, Kmove(knext))
            
            for ri in range(len(Rnext[0])):
                rnext = (Rnext[0][ri], Rnext[1][ri])
                if attacked[rnext]: # cannot move here
                    continue
                attacked = np.logical_and(attacked, Rmove(rnext))
                
                for bi in range(len(Bnext[0])):
                    bnext = (Bnext[0][bi], Bnext[1][bi])
                    if attacked[bnext]: # cannot move here
                        continue
                    attacked = np.logical_and(attacked, Bmove(bnext))
                    
                    for ni in range(len(Nnext[0])):
                        nnext = (Nnext[0][ni], Nnext[1][ni])
                        if attacked[nnext]: # cannot move here
                            continue
                        
                        moveIt(m, chessboard, Q, K, R, B, N,
                               qnext, knext, rnext, bnext, nnext)


# perform moves
def moveIt(m, chessboard, Q, K, R, B, N,
           qnext, knext, rnext, bnext, nnext):
    global largestpsum, optimalboard
    
    # push to stack
    Q.append(qnext)
    K.append(knext)
    R.append(rnext)
    B.append(bnext)
    N.append(nnext)
    
    # mark the board
    chessboard[Q[-1]] = m
    chessboard[K[-1]] = m
    chessboard[R[-1]] = m
    chessboard[B[-1]] = m
    chessboard[N[-1]] = m
    
    # start searching the next depth
    if m < 8:
        if m == 7: # check product after 7 moves
            rprods, cprods = boardProd(chessboard, Q, K, R, B, N)
            checkProd = np.logical_and(rprods == rprod, cprods == cprod).all()
        else:
            checkProd = True
        if checkProd:
            findMove(m + 1, chessboard, Q, K, R, B, N)
    
    # check final product sum
    if m == 8:
        rprods, cprods = boardProd(chessboard, Q, K, R, B, N)
        psum = sum(rprods) + sum(cprods)
        if psum > largestpsum:
            largestpsum = psum
            optimalboard = chessboard.copy()
            print(largestpsum)
            print(chessboard)
            print('')
    
    # pop from stack and unmark the board
    chessboard[Q.pop()] = -1
    chessboard[K.pop()] = -1
    chessboard[R.pop()] = -1
    chessboard[B.pop()] = -1
    chessboard[N.pop()] = -1


# calculate board products
def boardProd(chessboard, Q, K, R, B, N):
    # remove 0
    chessboard[Q[0]] = 1
    chessboard[K[0]] = 1
    chessboard[R[0]] = 1
    chessboard[B[0]] = 1
    chessboard[N[0]] = 1
    
    # -1 is taken care of by abs
    rprods = abs(np.prod(chessboard, axis=0))
    cprods = abs(np.prod(chessboard, axis=1))
    
    # add back 0
    chessboard[Q[0]] = 0
    chessboard[K[0]] = 0
    chessboard[R[0]] = 0
    chessboard[B[0]] = 0
    chessboard[N[0]] = 0
    
    return rprods, cprods


# initialization
n = 8 # board size
rprod = np.array([7, 1890, 8, 10080, 20, 840, 144, 1260]) # row products after 7 moves
cprod = np.array([2744, 36, 375, 336, 108, 240, 20, 504]) # column products after 7 moves
chessboard = np.ones((n,n))*-1 # board
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
with open('Puzzle_2016_Sep_result.txt','wb') as f:
    f.write(largestpsum)
    f.write('\n')
    f.write(str(optimalboard))