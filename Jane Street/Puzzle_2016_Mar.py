# approach the knigh moves puzzle with depth first search
import numpy, sys, copy

# define the knigh class and its operations
class Knight:
    def __init__(self):
        self.i = 0 # current row on the grid
        self.j = 0 # current column on the grid
        self.k = 0 # current step to mark the chessboard

    def jumpUL(self,forward): # forward = 1, jump forward, forward = -1, jump backward
        if forward != 1 and forward != -1:
            sys.stderr.write('Knight is not given correct time step')
            sys.exit(1)
        self.i -= 2
        self.j -= 1
        self.k += forward

    def jumpUR(self,forward):
        if forward != 1 and forward != -1:
            sys.stderr.write('Knight is not given correct time step')
            sys.exit(1)
        self.i -= 2
        self.j += 1
        self.k += forward

    def jumpLU(self,forward):
        if forward != 1 and forward != -1:
            sys.stderr.write('Knight is not given correct time step')
            sys.exit(1)
        self.i -= 1
        self.j -= 2
        self.k += forward

    def jumpRU(self,forward):
        if forward != 1 and forward != -1:
            sys.stderr.write('Knight is not given correct time step')
            sys.exit(1)
        self.i -= 1
        self.j += 2
        self.k += forward

    def jumpDL(self,forward):
        if forward != 1 and forward != -1:
            sys.stderr.write('Knight is not given correct time step')
            sys.exit(1)
        self.i += 2
        self.j -= 1
        self.k += forward

    def jumpDR(self,forward):
        if forward != 1 and forward != -1:
            sys.stderr.write('Knight is not given correct time step')
            sys.exit(1)
        self.i += 2
        self.j += 1
        self.k += forward

    def jumpLD(self,forward):
        if forward != 1 and forward != -1:
            sys.stderr.write('Knight is not given correct time step')
            sys.exit(1)
        self.i += 1
        self.j -= 2
        self.k += forward

    def jumpRD(self,forward):
        if forward != 1 and forward != -1:
            sys.stderr.write('Knight is not given correct time step')
            sys.exit(1)
        self.i += 1
        self.j += 2
        self.k += forward

def searchCurrent(knight): # check step validity
    global chessBoard
    # check if current position is valid
    if knight.i >= N or knight.j >= N or knight.i < 0 or knight.j < 0:
        sys.stderr.write('Knight fall off chess board')
        sys.exit(1)
    # check if current position is forbidden by preset marks or visited before
    currGrid = chessBoard[knight.i][knight.j]
    if currGrid != 0 and currGrid != knight.k:
##        print('forbidden at\n'+str(chessBoard))
        return
    # check if preset marks are violated
    if knight.k in presetMarks and currGrid != knight.k:
##        print('preset at\n'+str(chessBoard))
        return
    # check if the row and col sums are exceeded
    rs = 0
    for j in range(N):
        rs += chessBoard[knight.i][j]
    if currGrid == 0:
        rs += knight.k
    if rs > rowSums[knight.i]:
##        print('row sum violated at\n'+str(chessBoard))
        return
    cs = 0
    for i in range(N):
        cs += chessBoard[i][knight.j]
    if currGrid == 0:
        cs += knight.k
    if cs > colSums[knight.j]:
##        print('row sum violated at\n'+str(chessBoard))
        return
    # this step is valid
    searchNext(knight)

def searchNext(knight): # action
    global chessBoard, currPath
    # register current to path and mark the chess board
    currPath.append([knight.i, knight.j])
    currGrid = chessBoard[knight.i][knight.j]
    chessBoard[knight.i][knight.j] = knight.k
    # if finished, check score
    if knight.k == finishStep:
        checkScore(chessBoard, currPath)
    else:
    # knight jumps
        if knight.i > 1 and knight.j > 0:
            knight.jumpUL(1)
            searchCurrent(knight)
            knight.jumpDR(-1)
        if knight.i > 1 and knight.j < N-1:
            knight.jumpUR(1)
            searchCurrent(knight)
            knight.jumpDL(-1)
        if knight.i > 0 and knight.j > 1:
            knight.jumpLU(1)
            searchCurrent(knight)
            knight.jumpRD(-1)
        if knight.i > 0 and knight.j < N-2:
            knight.jumpRU(1)
            searchCurrent(knight)
            knight.jumpLD(-1)
        if knight.i < N-2 and knight.j > 0:
            knight.jumpDL(1)
            searchCurrent(knight)
            knight.jumpUR(-1)
        if knight.i < N-2 and knight.j < N-1:
            knight.jumpDR(1)
            searchCurrent(knight)
            knight.jumpUL(-1)
        if knight.i < N-1 and knight.j > 1:
            knight.jumpLD(1)
            searchCurrent(knight)
            knight.jumpRU(-1)
        if knight.i < N-1 and knight.j < N-2:
            knight.jumpRD(1)
            searchCurrent(knight)
            knight.jumpLU(-1)
    # remove changes
    currPath.pop()
    if currGrid == 0: # only remove marks not preset
        chessBoard[knight.i][knight.j] = 0
        
def checkScore(chessBoard, currPath):
    global maxScore, optPath
    print('reached checking stage at\n'+str(chessBoard))
    # first check spiral symmetry
    for i in range(N):
        for j in range(N-i*2-1):
            x = i
            y = i + j
            a = chessBoard[x][y]
            b = chessBoard[y][N-1-x]
            c = chessBoard[N-1-x][N-1-y]
            d = chessBoard[N-1-y][x]
            if not ((a == 0 and b == 0 and c == 0 and d == 0) or (a != 0  and b != 0 and c != 0 and d != 0)):
                return
    # then check rol and col sums
    for i in range(N):
        rs = 0
        for j in range(N):
            rs += chessBoard[i][j]
        if rs != rowSums[i]:
            return
    for j in range(N):
        cs = 0
        for i in range(N):
            cs += chessBoard[i][j]
        if cs != colSums[j]:
            return
    # finally check the largest product
    for i in range(N):
        currScore = 1
        for j in range(N):
            currGrid = chessBoard[i][j]
            if currGrid != 0:
                currScore = currScore * currGrid
        if currScore > maxScore:
            maxScore = currScore
            optPath = copy.deepcopy(currPath)
    for j in range(N):
        currScore = 1
        for i in range(N):
            currGrid = chessBoard[i][j]
            if currGrid != 0:
                currScore = currScore * currGrid
        if currScore > maxScore:
            maxScore = currScore
            optPath = copy.deepcopy(currPath)
    print(maxScore)
    print(str(optPath))

# initialize
print('Initializing ...')
N = 8 # grid size
finishStep = 28
chessBoard = numpy.zeros((N,N), dtype = int)
chessBoard[2][3] = 11
chessBoard[3][5] = 14
chessBoard[4][2] = 8
chessBoard[5][4] = 15
presetMarks = [8,11,14,15]
rowSums = [10,34,108,67,63,84,24,16]
colSums = [7,14,72,66,102,90,42,13]
maxScore = 1 # maximum score found so far
currPath = [] # current path in the form: [[i,j],...]
optPath = [] # path corresponding to maximum score
knight = Knight()

for i in range(N):
    for j in range(N):
        print('Searching starting position: '+str(i)+','+str(j))
        knight.i = i
        knight.j = j
        knight.k = 1
        searchCurrent(knight)

# save the result
import pickle
f = open('Puzzle_2016_Mar_result.pckl','w')
pickle.dump(maxScore,f)
pickle.dump(optPath,f)
f.close()

f = open('Puzzle_2016_Mar_result.txt','w')
f.write('Score: '+str(maxScore)+'\n\n')
f.write('Path:\n')
# produce the optimal chess board
optCB = numpy.zeros((N,N), dtype = int)
for i in range(len(optPath)):
    optCB[optPath[i][0]][optPath[i][1]] = i+1
f.write(str(optCB))
f.close()
